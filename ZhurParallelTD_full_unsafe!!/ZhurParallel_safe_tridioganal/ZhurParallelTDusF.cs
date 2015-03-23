using System;
using System.Threading;
using System.Threading.Tasks;
using System.Text;
using System.Collections.Concurrent;
using System.Runtime.InteropServices;

namespace ZhurParallelTDusF
{
    public unsafe class EPlTDusF
    {
        
        #region EigenValueAndVector2x2
        /// <summary>
        /// Метод, рассчитывающий собственные числа и вектора матрицы 2х2.
        /// </summary>
        /// <param name="a">Указатель на исходную матрицу.</param>
        /// <param name="lambda">Указатель на массив собственных чисел по возрастанию.</param>
        /// <param name="vector">Указатель на матрицу собственных векторов.</param>
        static void EigenValueAndVector2x2(double* a, double* lambda, double* vector)
        {
            //получаем собственные значения
            *lambda = (*a + a[3]) * 0.5 - Math.Sqrt(((*a - a[3]) * (*a - a[3])) * 0.25 + Sqr(a[1]));
            lambda[1] = (*a + a[3]) * 0.5 + Math.Sqrt(((*a - a[3]) * (*a - a[3])) * 0.25 + Sqr(a[1]));

            //Создаем матрицу собственных векторов
            *vector = a[1];
            vector[1] = a[1];
            vector[2] = *lambda - *a;
            vector[3] = lambda[1] - *a;

            //Нормируем ее
            NormVectorMatrix2x2(vector);
        }

        /// <summary>
        /// Метод, рассчитывающий собственные числа и вектора матрицы 2х2, в случае трехдиагональной матрицы.
        /// </summary>
        /// <param name="a">Указатель на исходную трехдиагональную матрицу.</param>
        /// <param name="lambda">Указатель на массив собственных чисел по возрастанию.</param>
        /// <param name="vector">Указатель на последнюю строку матрицы собственных векторов.</param>
        static void EigenValueAndVector2x2TD(double* a, double* lambda, double* vector, bool up)
        {
            //получаем собственные значения
            *lambda = (*a + a[3]) * 0.5 - Math.Sqrt(((*a - a[3]) * (*a - a[3])) * 0.25 + Sqr(a[1]));
            lambda[1] = (*a + a[3]) * 0.5 + Math.Sqrt(((*a - a[3]) * (*a - a[3])) * 0.25 + Sqr(a[1]));

            //Создаем матрицу собственных векторов
            double[,] vtemp = new double[2, 2];
            fixed (double* v = vtemp)
            {
                *v = a[1]; 
		        v[1] = a[1];
                v[2] = *lambda - *a;
                v[3] = lambda[1] - *a;

                //Нормируем ее
                NormVectorMatrix2x2(v);

                //Получаем первую или последнюю строку
                if (up)
                {
                    *vector = v[0];
                    vector[1] = v[1];
                }
                else
                {
                    *vector = v[2];
                    vector[1] = v[3];
                }
            }
        }
        #endregion

        #region NormVectorMatrix2x2
        /// <summary>
        /// Метод, производящий нормировку матрицы собственных векторв.
        /// </summary>
        /// <param name="vector">Указатель на матрицу собственных векторов.</param>
        /// <param name="n">Размерность матрицы.</param>>
        static void NormVectorMatrix2x2(double* vector)
        {
            Parallel.For(0, 2, j =>
            {
                double temp = 0.0d;

                for (double* start = &vector[j], stop = start + 4; start < stop; start += 2)
                    temp = temp + Sqr(*start);

                temp = Math.Sqrt(temp);

                for (double* start = &vector[j], stop = start + 4; start < stop; start += 2)
                {
                    /***********************************************
                    * В зависимости от используемого процессора можно 
                    * выбрать каким способом считать обратный корень.
                    * Если Процессор имеет большее количество модулей
                    * по работе с плавающей запятой, чем с целыми числами,
                    * стоит выбрать стандартный метод Math.Sqrt(x),
                    * в ином случае метод InvSqrt(x) будет работать
                    * быстрее, хоть и с небольшой погрешностью.
                    * При равенстве модулей счета - производительность 
                    * первого или второго метода нужно проверять отдельно. 
                    *************************************************/
                    *start = *start / temp; // * InvSqrt(temp);                        
                }
            });
        }

        #endregion

        #region NormVectorMatrixBD
        /// <summary>
        /// Нормирует матрицу векторов.
        /// </summary>
        /// <param name="d">Указатель на матрицу векторов.</param>
        /// <param name="n">Размерность матрицы векторов.</param>
        static void NormVectorMatrixBD(BigDouble* d, int n)
        {
            BigDouble[,] res = new BigDouble[n, n];

            Parallel.For(0, n, j =>
            {
                BigDouble temp = new BigDouble(0.0, 1);

                for (BigDouble* st = d + j, sp = st + (n - 1) * n + 1; st < sp; st += n)
                    temp += BigDouble.Pow(*st, 2);

                temp = BigDouble.Sqrt(temp);

                for (BigDouble* st = d + j, sp = st + (n - 1) * n + 1; st < sp; st += n)
                    *st = *st / temp;
            });

        }
        #endregion

        #region FormANk
        /// <summary>
        /// Формируем матрицу ANk исходя из ANk = Inverse(Unk) * ANk * UNk.
        /// </summary>
        /// <param name="v">Указатель на последнюю строку матрицы собственных векторов.</param>
        /// <param name="d">Указатель на главную диагональ исходной матрицы.</param>
        /// <param name="e">Указатель на побочную диагональ исходной матрицы.</param>
        /// <param name="a">Указатель на ненулевой столбец ограниченно-диагональной матрицы.</param>
        /// <param name="l">Указатель на массив собственных чисел.</param>
        /// <param name="n">Размерность матрицы.</param>
        static void FormANk(double* v, double* d, double* e, double* a, double* l, int n)
        {
            //Получаем множитель ненулевого столбца
            double ak = e[n - 1];

            //Зануляем (n-3)й элемент побочной диагонали
            e[n - 2] = 0.0d;

            //Основной цикл формирования ограниченно-диагональной матрицы
            Parallel.For(0, n, i =>
            {
                //Формируем главную диагональ
                d[i] = l[i];
                //Формируем ненулевой столбец
                a[i] = ak * v[i];
            });
        }

        /// <summary>
        /// Формирует ограниченно-диагональную матрицу ANk для последнего шага.
        /// </summary>
        /// <param name="AN">Указатель на массив, составляющий последний (новый) столбец исходной матрицы.</param>
        /// <param name="V">Указатель на матрицу собственных векторов матрицы размерности N-1, вырезанной из исходной.</param>
        /// <param name="tempA">Указатель на временный массив.</param>
        /// <param name="N">Размерность матрицы.</param>
        static void FormANkAdd(double* AN, double* V, double* tempA, int N)
        {
            Parallel.For(0, N, i =>
            {
                var temp = 0.0d;
                for (int m = 0; m < N; m++)
                    temp += AN[m] * V[m * N + i];
                tempA[i] = temp;
            });

            tempA[N] = AN[N];
        }
        #endregion

        #region EigenValue
        /// <summary>
        /// Метод рассчитывает собственные вектора и числа заданной ограниченнно-диагональной матрицы с указанной точностью.
        /// </summary>
        /// <param name="D">Указатель на массив, составляющий главную диагональ исходной матрицы.</param>
        /// <param name="E">Указатель на массив, составляющий побочную диагональ исходной матрицы.</param>
        /// <param name="A">Указатель на массив, составляющий ненулевой столбец ограниченно-диагональной матрицы.</param>
        /// <param name="n">Размерность ограниченно-диагональной матрицы (или шаг вычисления).</param>
        /// <param name="N">Размерность исходной матрицы</param>
        /// <param name="lambda">Указатель на массив, составляющий сосбтвенные числа.</param>
        /// <param name="vector">Указатель на массив, составляющий последнюю строку матрицы собственных векторов.</param>
        /// <param name="epsilon">Точность вычисления.</param>
        static void EigenValue(double* D, double* E, double* A, int n, int N, double* lambda, double* vector, double epsilon)
        {

            //Находим собственные значения на внутренних отрезках.                        
            Parallel.For(0, n - 2, i =>            
            {
                lambda[i + 1] = EigenValueInterval(D, A, i, n, epsilon);
            });

            //Находим собственные значения на граничных интервалах.     
            EigenValueInfinit(D, A, lambda, n - 1);            
            
            //Находим собственные вектора исходной
            //ограниченно-диагональной матрицы.            
            EigenVector(vector, lambda, A, D, n);
        }

        /// <summary>
        /// Метод рассчитывает собственные вектора и числа заданной ограниченнно-диагональной матрицы с указанной точностью.
        /// </summary>
        /// <param name="D">Указатель на массив, составляющий главную диагональ исходной матрицы.</param>
        /// <param name="E">Указатель на массив, составляющий побочную диагональ исходной матрицы.</param>
        /// <param name="A">Указатель на массив, составляющий ненулевой столбец ограниченно-диагональной матрицы.</param>
        /// <param name="n">Размерность ограниченно-диагональной матрицы (или шаг вычисления).</param>
        /// <param name="N">Размерность исходной матрицы</param>
        /// <param name="lambda">Указатель на массив, составляющий сосбтвенные числа.</param>
        /// <param name="vector">Указатель на массив, составляющий последнюю строку матрицы собственных векторов.</param>
        /// <param name="epsilon">Точность вычисления.</param>
        static void EigenValueAdd(double* D, double* E, double* A, int N, double* lambda, double* vector, double epsilon)
        {
            //Находим собственные значения на внутренних отрезках.                        
            Parallel.For(0, N - 2, i =>
            {
                lambda[i + 1] = EigenValueInterval(D, A, i, N, epsilon);
            });

            //Находим собственные значения на граничных интервалах.                        
            EigenValueInfinit(D, A, lambda, N - 1);

            //Находим собственные вектора исходной
            //ограниченно-диагональной матрицы.
            EigenVectorAdd(vector, lambda, A, D, N);
        }
        #endregion

        #region EigenValueInfinit
        /// <summary>
        /// Находит краевые собственные значения ограниченнно-диагональной матрицы.
        /// </summary>
        /// <param name="d">Указатель на массив диагональных элементов ограниченнно-диагональной матрицы.</param>
        /// <param name="A">Указатель на ненулевой столбец ограниченнно-диагональной матрицы.</param>
        /// <param name="Lambda">Указатель на массив сосбтвенных чисел ограниченнно-диагональной матрицы.</param>
        /// <param name="n">Шаг вычисления.</param>
        static void EigenValueInfinit(double* d, double* A, double* Lambda, int n)
        {
            double D = d[n];

            //SA
            double SA = 0.0d;            
            for (double* i = A, i_end = A + n; i < i_end; i++)
                SA += Sqr(*i);

            //Sl
            double Sl = 0.0d;            
            for (double* i = d, i_end = d + n; i < i_end; i++)
                Sl += *i;

            //SL
            double SL = 0.0d;
            for (double* i = Lambda + 1, i_end = Lambda + n; i < i_end; i++)
                SL += *i;

            //Sli
            double Sli = 0.0d;
            for (double* i = d, i_end = d + n - 1; i < i_end; i++)
            {
                double temp = 0.0d;
                for (double* j = i + 1, j_end = d + n; j < j_end; j++)
                    temp += *j;
                Sli += *i * temp;
            }

            //SLi
            double SLi = 0.0d;
            for (double* i = Lambda, i_end = Lambda + n - 1; i < i_end; i++)
            {
                double temp = 0.0d;
                for (double* j = i + 1, j_end = Lambda + n; j < j_end; j++)
                    temp += *j;
                SLi += *i * temp;
            }

            //S
            double S = D + Sl - SL;

            //Дескриминант квадратного уравнения
            double squareroot = Math.Sqrt(S * S + 4 * (SL * S + SLi - D * Sl - Sli + SA));

            *Lambda = (S - squareroot) * 0.5d;
            Lambda[n] = (S + squareroot) * 0.5d;
        }

        #endregion

        #region EigenVectors
        /// <summary>
        /// Находит собственные вектора симметричной трехдиагональной матрицы по известным собственным значениям.
        /// </summary>
        /// <param name="D">Указатель на массив, составляющий главную диагональ трехдиагональной матрицы.</param>
        /// <param name="E">Указатель на массив, составляющий побочную диагональ трехдиагональной матрицы.</param>
        /// <param name="L">Указатель на массив собственных чисел.</param>
        /// <param name="n">Размерность матрицы.</param>
        /// <param name="res">Указатель на матрицу собственных векторов трехдиагональной матрицы.</param>
        private static void EigenVectors(double* D, double* E, double* L, int n, double* res)
        {
            //По скольку, при вычислении собственных векторов
            //трехдиагональной матрицы, значения временно могут
            //выходить за рамки значений типа Double,
            //создаем матрицу специального типа.
            BigDouble[,] vector = new BigDouble[n, n];

            //Вычисляем собственные вектора
            //трехдиагональной матрицы в BigDouble
            Parallel.For(0, n, j =>
            {
                vector[0, j] = new BigDouble(1.0d, -100);
                vector[1, j] = ((BigDouble)(L[j] - *D)) * vector[0, j] / (BigDouble)(*E);


                double* d = D + 1, e = E;
                fixed (BigDouble* v = &vector[2, j])
                    for (BigDouble* st = v, sp = v + (n - 3) * n + 1; st < sp; st += n, d++, e++)
                    {
                        *st = ((BigDouble)(L[j] - *d) * *(st - n) - (BigDouble)(*e) * *(st - n - n)) / (BigDouble)(*(e + 1));
                    }
            });

            //Нормируем вектора
            fixed (BigDouble* v = &vector[0, 0])
                NormVectorMatrixBD(v, n);

            //Преобразуем полученные значения
            //к станартному типу Double и 
            //заполняем соответствующую матрицу.
            Parallel.For(0, n, i =>
            {
                double* r = res + i * n;
                fixed (BigDouble* v = &vector[i, 0])
                    for (BigDouble* st = v, sp = st + n; st < sp; st++, r++)
                    {
                        *r = (double)(*st);
                    }
            });
        }
        #endregion        

        #region EigenVector
        /// <summary>
        /// Метод вычисляет собственные вектора исходной ограниченно-диагональной матрицы.
        /// </summary>
        /// <param name="V">Указатель на последнюю строку матрицы собственных векторов ограниченно-диагональной матрицы.</param>
        /// <param name="lambda">Указатель на массив собственных чисел ограниченно-диагональной матрицы.</param>
        /// <param name="A">Указатель на ненулевой столбец ограниченно-диагональной матрицы.</param>
        /// <param name="D">Указатель на массив, составляющий главную диагональ матрицы.</param>
        /// <param name="n">Размерность матрицы</param>
        static void EigenVector(double* V, double* lambda, double* A, double* D, int n)
        {
            //var bags = new ConcurrentBag<double[]>();
            //Вычисляем собственные вектора и
            //заносим последнюю строку в массив.
            Parallel.For(0, n, j =>            
            {
                //V[j] = eigenVectorSubroutine(lambda[j], A, D, n, lambda);
                double temp = 1.0d;

                for (double* i = D, _a = A, i_end = D + n - 1; i < i_end; i++, _a++)
                {
                    var t = Sqr(*_a / (lambda[j] - *i));
                    temp += t;
                }

                V[j] = Math.Sqrt(1.0d / temp);
                
            });
            
        }

        private static double eigenVectorSubroutine(double lambda, double* A, double* D, int n, double* l)
        {
            double* vector = stackalloc double[n - 1];
            
            var _i = -1;

            for (int i = 0; i < n - 1; i++)
                if (lambda != D[i])
                    vector[i] = A[i] / (lambda - D[i]);
                else
                    _i = i;
            
            if (_i != -1)
            {
                var temp = D[n - 1] - lambda;

                for (int i = 0; i < n - 1; i++)
                    if (i != _i)
                        temp += A[i] * vector[i];

                vector[_i] = temp * (-1.0d / A[_i]);
            }
            
            var t = 1.00d;

            for (double* i = vector, i_end = vector + n - 1; i < i_end; i++)
                t += Sqr(*i);

            return Math.Sqrt(1.0d / t); 
        }

        /// <summary>
        /// Метод вычисляет собственные вектора исходной ограниченно-диагональной матрицы.
        /// </summary>
        /// <param name="V">Указатель на последнюю строку матрицы собственных векторов ограниченно-диагональной матрицы.</param>
        /// <param name="lambda">Указатель на массив собственных чисел ограниченно-диагональной матрицы.</param>
        /// <param name="A">Указатель на ненулевой столбец ограниченно-диагональной матрицы.</param>
        /// <param name="D">Указатель на массив, составляющий главную диагональ матрицы.</param>
        /// <param name="n">Размерность матрицы</param>
        static void EigenVectorAdd(double* V, double* lambda, double* A, double* D, int n)
        {
            //Вычисляем собственные вектора и
            //заносим последнюю строку в массив.
            Parallel.For(0, n, j =>
            {
                var temp = 1.0d;
                
                *(V + j + (n - 1) * n) = 1.0d;

                for (double* i = D, _a = A, i_end = D + n - 1, v_st = V + j; i < i_end; i++, _a++, v_st += n)
                {
                    *v_st = *_a / (lambda[j] - *i);
                    temp += Sqr(*v_st);
                }

                for (double* v_st = V + j, v_end = V + j + n * (n - 1) + 1; v_st < v_end; v_st += n)
                {
                    *v_st = Math.Sqrt(*v_st / temp);                    
                }                

            });

        }
        #endregion

        #region DiagonalVector
        /// <summary>
        /// Внутренняя процедура вычисления собственных векторов.
        /// </summary>
        /// <param name="vector">Указатель на матрицу векторов матрицы N-1.</param>
        /// <param name="tmp">Указатель на временную матрицу.</param>
        /// <param name="N">Размерность матрицы.</param>
        static void DiagonalVector(double* vector, double* tmp, int N)
        {
            Parallel.For(0, N - 1, i =>
            {
                for (int j = 0; j < N - 1; j++)
                    tmp[i * N + j] = vector[i * (N - 1) + j];
            });

            tmp[N * N - 1] = 1.0d;
        }
        #endregion

        #region MatrixMultiply
        /// <summary>
        /// Метод перемножения двух квадратных матриц.
        /// </summary>
        /// <param name="M1">Указатель на первую матрицу.</param>
        /// <param name="M2">Указатель на вторую матрицу.</param>
        /// <param name="result">Указатель на полученную матрицу.</param>
        /// <param name="n">Размерность матриц.</param>
        static void MatrixMultiply(double* M1, double* M2, double* result, int n)
        {
            Parallel.For(0, n, i =>
            {
                for (double* k_start = M1 + i * n, k_end = k_start + n, j_end = result + i * n + n, y_tmp = M2; k_start < k_end; k_start++)
                {
                    //сохраняем текущее значение для многократного использования
                    double temp = *k_start;
                    for (double* _ij = result + i * n; _ij < j_end; _ij++, y_tmp++)
                        *_ij += temp * *y_tmp;
                }
            });
        }    
        #endregion

        #region vectorsAdd
        /// <summary>
        /// Заполняем массивы главной диагонали, составляющей собственные числа, и
        /// последнего (нового) столбца матрицы.
        /// </summary>
        /// <param name="AN">Указатель на массив, составляющий последний (новый) столбец матрицы.</param>
        /// <param name="A">Указатель на исходную матрицу.</param>
        /// <param name="lambda">Указатель на массив собственных чисел.</param>
        /// <param name="D">Указатель на массив, составляющий главную диагональ.</param>
        /// <param name="N">Размерность матрицы.</param>
        static void vectorsAdd(double* AN, double* A, double* lambda, double* D, int N)
        {
            //Заполняем массив последнего столбца
            Parallel.For(0, N, i =>
            {
                AN[i] = A[i * N + N - 1];
            });

            //Заполняем главную диагональ.
            Parallel.For(0, N - 1, i =>
            {
                D[i] = lambda[i];
            });
        }
        #endregion


        #region EigenValueInterval
        /// <summary>
        /// Метод вычисляет собственное значение на заданном интервале.
        /// </summary>
        /// <param name="d">Указатель на массив, составляющий главную диагональ матрицы.</param>
        /// <param name="a">Указатель на массив, составляющий ненулевой столбец матрицы.</param>
        /// <param name="k">Номер отрезка.</param>
        /// <param name="n">Размерность матрицы.</param>
        /// <param name="epsilon">Точность вычисления.</param>
        /// <returns>Найденное собственное значение.</returns>
        static double EigenValueInterval(double* d, double* a, int k, int n, double epsilon)
        {
            //Получаем значение
            //точек интервала.
            double L1 = d[k];
            double L2 = d[k + 1];

            double lambda = 0.0d;

            var tempA = stackalloc double[n - 1];
            for (double* start = a, stop = a + n - 1, _a = tempA; start < stop; start++, _a++)
                *_a = Sqr(*start);

            //Получаем максимальное значение делений на 2,
            //достаточных для нахождения собственного числа
            //с заданной точностью:
            // N = log([L2 - L1] / epsilon) / log(2)
            double _N = Math.Log((L2 - L1) / epsilon);
            // 1 / log(2) ~ 3.32192809488736
            int N = (int)Math.Ceiling(_N * 3.32192809488736d);


            //Если заданная точность больше значения
            //точек интервала (N < 0)
            if (N < 0)
            {
                for (int x = 0; x < (int)Math.Ceiling(-_N) + 3; x++)
                    epsilon *= 0.1d;

                //Основной цикл нахождения собственного 
                //значения на заданном интервале.
                for (int j = 0; j < -N + 1; j++)
                {
                    //Предполагаемое значение.
                    lambda = (L1 + L2) * 0.5d;

                    double temp = d[n - 1];

                    //Подставляем значение lambda в уравнение G(λ) = D + (A1*A1)/(λ-λ1) + (A2*A2)/(λ-λ2) +...+ (Ak*Ak)/(λ-λk),
                    //где λi - собственное значение матрицы Bk, Ai - значение массива A.  

                    for (double* _a = tempA, _a_end = tempA + n - 1, _d = d; _a < _a_end; _a++, _d++)
                        temp += *_a / (lambda - *_d);

                    //Если |G(λ) - λ| < epsilon - собственное число с точностью epsilon
                    //найдено - завершение цикла.
                    if (Math.Abs(temp - lambda) <= epsilon)
                    {
                        break;
                    }

                    //Иначе уменьшаем интервал в двое вправо или влево.
                    if (temp - lambda > 0.0d)
                        L1 = lambda;
                    else
                        L2 = lambda;
                }
            }
            else
                //Основной цикл нахождения собственного 
                //значения на заданном интервале.
                for (int j = 0; j < N + 1; j++)
                {
                    //Предполагаемое значение.
                    lambda = (L2 + L1) * 0.5d;

                    double temp = d[n - 1];

                    //Подставляем значение lambda в уравнение G(λ) = D + (A1*A1)/(λ-λ1) + (A2*A2)/(λ-λ2) +...+ (Ak*Ak)/(λ-λk),
                    //где λi - собственное значение матрицы Bk, Ai - значение массива A.    

                    for (double* _a = tempA, _a_end = tempA + n - 1, _d = d; _a < _a_end; _a++, _d++)
                        temp += *_a / (lambda - *_d);

                    //Если |G(λ) - λ| < epsilon - собственное число с точностью epsilon
                    //найдено - завершение цикла.
                    if (Math.Abs(temp - lambda) <= epsilon)
                        break;

                    //Иначе уменьшаем интервал в двое вправо или влево.
                    if (temp - lambda > 0.0d)
                        L1 = lambda;
                    else
                        L2 = lambda;
                }

            //Возвращаем найденное собственное значение.
            return lambda;
        }
        #endregion

        #region Sqr
        /// <summary>
        /// Возвращает квадрат заданного числа.
        /// </summary>
        /// <param name="p">Исходное число.</param>
        /// <returns>Квадрат числа.</returns>
        static double Sqr(double p)
        {
            return p * p;
        }
        #endregion

        #region Fast Inerse Square Root
        /// <summary>
        /// Метод быстрого обратного корня. 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
      /*  static unsafe double InvSqrt(double x)
        {
            double xhalf = 0.5 * x;
            long i = *(long*)&x;
            i = 0x5fe6ec85e7de30da - (i >> 1); //5f375a86 для int и float
            x = *(double*)&i;
            x = x * (1.5f - xhalf * x * x);
            x = x * (1.5f - xhalf * x * x);
            x = x * (1.5f - xhalf * x * x); //Повторяя данную строку точность вычисления увеличивается.
            return x;
        }*/
        #endregion


        #region smatrixevv
        /// <summary>
        /// Метод вычисляет собственные вектора и собственные значения исходной симметричной матрицы A 
        /// c заданной точностью epsilon, путем преобразования ее к трехдиагональному виду.
        /// Алгоритм основан на <see cref="">Итерационном методе вычисления собственных значений симметричных матриц авторства Журавлев В. М., Журавлев А. В.</see>
        /// </summary>
        /// <param name="A">Исходная матрица.</param>
        /// <param name="lambda">Собственные числа.</param>
        /// <param name="vector">Собственные вектора.</param>
        /// <param name="epsilon">Точность вычисления.</param>
        /// <param name="vectorsFound">Индикатор, указывающий на надобность вычисления собственных векторов.
        /// <value>false</value> - если вычисление собственных векторов не требуется и <value>true</value> - в ином случае.</param>
        /// <returns>True, если собственные вектора и числа найдены, иначе - false.</returns>       
        public static bool smatrixevv(double[,] A, ref double[] lambda, ref double[,] vector, double epsilon, bool vectorsFound)
        {
            //Создаем начальные массив собственных чисел
            //и матрицу собственных векторов.
            lambda = new double[2];
            vector = new double[2, 2];

            //Получаем размерность исходной матрицы.
            int n = A.GetLength(0);

            //Если исходная матрица размерности 2х2,
            //то вычисляем собственные значения и вектора 
            //и выходим из метода.            
            if (n == 2)
            {
                fixed (double* l = lambda, v = vector, _a = A)
                EigenValueAndVector2x2(_a, l, v);
                return true;
            }
            
            //Подготовка к преобразованию матрицы
            //к трехдиагональному виду. Создаем
            //массивы главной, побочной диагоналей
            //и факторов H(i), а так же матрицу
            //преобразования.
            double[] dA = new double[n];
            double[] eA = new double[n - 1];
            double[] tau = new double[1];
            double[,] Q;

            //Преобразование симметричной матрицы к трехдиагональному
            //виду. Используется метод из свободной библиотеки AlgLib.
            smatrix2TD.smatrixtd(ref A, n, true, ref tau, ref dA, ref eA);

            //Нахождаение собственных чисел и векторов трехдиагональной матрицы.
            bool result = smatrixevvTD(dA, eA, n, ref lambda, ref vector, epsilon, vectorsFound);
            
            //Получаем собственные вектора
            //исходной матрицы А, если требуется.
            if (vectorsFound)
            {                
                Q = new double[n, n];
                //Распаковываем матрицу преобразования.
                //Используется метод из свободной библиотеки AlgLib.
                smatrix2TD.smatrixtdunpackq(A, n, true, tau, ref Q); 

                double[,] temp = (double[,])vector.Clone();
                vector = new double[n, n];
                //Вычисляем собственные вектора исходной матрицы
                //путем произведения матрицы преобрзаования на 
                //матрицу собствнных векторов трехдиагональной матрицы.
                fixed (double* q = Q, v = vector, t = temp)
                    MatrixMultiply(q, t, v, n);

            }

            return result;
        }
        #endregion

        public static bool smatrixTDtemp(double[] dA, double[] eA, int n, ref double[] lambda,
                                        ref double[,] vector, double epsilon, bool vectorsFound)
        {
            double[] d = (double[])dA.Clone();
            double[] e = (double[])eA.Clone();

            //Текущий размер вычисляемых матриц
            int del = 2;
            //Количество матриц текущего размера
            int ni = (n + 1) % (del + 1);
            //Количество оставшихся элементов
            int tale = n + 1 - ni * (del + 1);
            //Массив верхних или нижних строк собственых
            //векторов матриц текущего размера
            double[][] vectors = new double[ni][];
            double[][] lambdas = new double[ni][];

            //Вычисляем собственные числа и вектора
            //всех матриц, размерности 2
            for(int i = 0, j = 1; i < n; i += (del + 1), j++)
            {
                double[,] Ak = new double[del, del];
                Ak[0, 0] = dA[i];
                Ak[0, 1] = eA[i];
                Ak[1, 0] = eA[i];
                Ak[1, 1] = dA[i + 1];

                double[] L = new double[del];
                double[] V = new double[del];
                fixed (double* a = Ak, l = L, v = V)
                    EigenValueAndVector2x2TD(a, l, v, j % 2 == 0);

                lambdas[j] = L;
                vectors[j] = V;
            }

            //Если количество матриц не четное,
            //то нужно совместить последнюю матрицу
            //размерности del с остатком
            if (ni % 2 != 0 && tale > 0)
            {                
                //Посчитать остатток методом Add
                for (int i = 1; i <= tale; i++)
                {
                    double[] L = new double[del + i];
                    double[] V = new double[del + i];
                    fixed (double* l = L, v = V, l_old = lambdas[ni - 1], v_old = vectors[ni - 1])
                        EigenValueAndVectorAddTemp(l, l_old, v, v_old, e[n - tale - 2 + i], 
                                                   d[n - tale - 1 + i], epsilon);
                }
            }

            //Обновляем значения
            del = del + del + 1;
            ni = (n + 1) % (del + 1);
            tale = n + 1 - ni * (del + 1);

            while (del <= n)
            {
                double[][] newVectors = new double[ni][];
                double[][] newLambdas = new double[ni][];

                int k = 0;

                for (int i = 0, j = 0; i < n; i += (del + 1), j++, k += 2)
                {
                    //Вычисление чисел и векторов
                    //Если j % 2 = 0, записываем первую строку
                    //матрицы сосбтвенных векторов, иначе - последнюю
                    int poz = del / 2 + i;
                    newLambdas[j] = new double[del];
                    newVectors[j] = new double[del];
                    fixed (double* oldL1 = lambdas[k], oldV1 = vectors[k],
                           oldL2 = lambdas[k + 1], oldV2 = vectors[k + 1],
                           newL = newLambdas[j], newV = newVectors[j])
                        EigenValueAndVectorXxX(d[poz], e[poz - 1], i, del, oldL1, oldV1, oldL2, 
                                               oldV2, newL, newV, epsilon, (j + 1) % 2 == 0);                   

                }

                //Если количество матриц не четное,
                //то нужно совместить последнюю матрицу
                //размерности del с остатком
                if (ni % 2 != 0 && tale > 0)
                    if (tale < 3)
                    {
                        //Посчитать остатток методом Add
                        for (int i = 1; i <= tale; i++)
                        {
                            double[] L = new double[del + i];
                            double[] V = new double[del + i];
                            fixed (double* l = L, v = V, l_old = lambdas[ni - 1], v_old = vectors[ni - 1])
                                EigenValueAndVectorAddTemp(l, l_old, v, v_old, e[n - tale - 2 + i],
                                                           d[n - tale - 1 + i], epsilon);
                        }
                    }
                    else
                    {
                        double[] L = new double[del + tale];
                        double[] V = new double[del + tale];
                        //Посчитать как 
                        //проверить значение k в этом месте!!!!
                        fixed (double* oldL = lambdas[k], oldV = vectors[k],
                            newL = newLambdas[ni - 1], newV = newVectors[ni - 1])
                            EigenValueAndVectorXxY(d[n - tale], e[n - tale - 1], del, tale - 1,
                                                   oldL, oldV, newL, newV, epsilon);
                    }

                lambdas = newLambdas;
                vectors = newVectors;

                //Обновляем значения
                del = del + del + 1;
                ni = (n + 1) % (del + 1);
                tale = n + 1 - ni * (del + 1);
            }

            lambda = lambdas[0];

            if (vectorsFound)
            {
                //Создаем матрицу собственных векторов.
                vector = new double[n, n];
                //Вычисляем собственные вектора
                //трехдиагональной матрицы.
                fixed (double* _d = dA, _e = eA, l = lambda, v = vector)
                    EigenVectors(_d, _e, l, n, v);
            }
            
            return true;
        }

        #region smatrixevvTD
        /// <summary>
        /// Метод вычисляет собственные числа и вектора трехдиагональной матрицы, 
        /// представленной массивами, составляющими главную и побочную диагональ трехдиагональной матрицы.
        /// </summary>
        /// <param name="dA">Массив, составляющий главную диагональ трехдиагональной матрицы.</param>
        /// <param name="eA">Массив, составляющий побочную диагональ трехдиагональной матрицы.</param>
        /// <param name="n">Размерность матрицы.</param>
        /// <param name="lambda">Собственные числа трехдиагональной матрицы.</param>
        /// <param name="vector">Собственные вектора трехдиагональной матрицы.</param>
        /// <param name="epsilon">Точность вычисления.</param>
        /// <param name="vectorsFound">Индикатор, указывающий на надобность вычисления собственных векторов.
        /// <value>false</value> - если вычисление собственных векторов не требуется и <value>true</value> - в ином случае.</param>
        /// <returns>True, если собственные вектора и числа найдены, иначе - false.</returns>
        public static bool smatrixevvTD(double[] dA, double[] eA, int n, ref double[] lambda,
                                        ref double[,] vector, double epsilon, bool vectorsFound)
        {
            //Создаем матрицу Ak, размерностью 2х2.
            double[,] Ak = new double[2, 2];
            Ak[0, 0] = dA[0];
            Ak[0, 1] = eA[0];      
            Ak[1, 0] = eA[0];
            Ak[1, 1] = dA[1];

            //Создаем переменные для хранения
            //промежуточных значений.
            double[] L = new double[2];
            double[] V = new double[2];

            //Находим собственные вектора и 
            //значения матрицы Ak.
            fixed (double* _a = Ak, l = L, v = V)
                EigenValueAndVector2x2TD(_a, l, v, false);

            //Создаем матрицу ANk.
            double[] d = (double[])dA.Clone();
            double[] e = (double[])eA.Clone();
            double[] a;

            //Основной цикл нахождения собственных
            //векторов и значений исходной матрицы.            
            for (int i = 2; i < n; i++)
            {
                //Создаем актуальную матрицу ANk.
                a = new double[i];
                fixed (double* _v = V, _d = d, _e = e, _a = a, l = L)
                    FormANk(_v, _d, _e, _a, l, i); //ANk = Inverse(UNk) * ANk * UNk;

                //Находим собственные вектора и значения 
                //матрицы Bk+1.
                L = new double[i + 1];
                V = new double[i + 1];
                fixed (double* l = L, v = V, _e = e, _a = a, _d = d)
                    EigenValue(_d, _e, _a, i + 1, n, l, v, epsilon);
            }            

            //Получаем собственные вектора
            //трехдиагональной матрицы, если требуется.
            if (vectorsFound)
            {
                //Создаем матрицу собственных векторов.
                vector = new double[n, n];
                //Вычисляем собственные вектора
                //трехдиагональной матрицы.
                fixed (double* _d = dA, _e = eA, l = L, v = vector)
                    EigenVectors(_d, _e, l, n, v);
            }

            lambda = L;

            return true;
        }

        /// <summary>
        /// Метод вычисляет собственные числа и вектора трехдиагональной матрицы, 
        /// представленной двумерным массивом.
        /// </summary>
        /// <param name="dA">Массив, составляющий главную диагональ трехдиагональной матрицы.</param>
        /// <param name="eA">Массив, составляющий побочную диагональ трехдиагональной матрицы.</param>        
        /// <param name="lambda">Собственные числа трехдиагональной матрицы.</param>
        /// <param name="vector">Собственные вектора трехдиагональной матрицы.</param>
        /// <param name="epsilon">Точность вычисления.</param>
        /// <param name="vectorsFound">Индикатор, указывающий на надобность вычисления собственных векторов.
        /// <value>false</value> - если вычисление собственных векторов не требуется и <value>true</value> - в ином случае.</param>
        /// <returns>True, если собственные вектора и числа найдены, иначе - false.</returns>
        public static bool smatrixevvTD(double[,] A, ref double[] lambda,
                                        ref double[,] vector, double epsilon, bool vectorsFound)
        {
            //Получаем размерность матрицы.
            int n = A.GetLength(0);

            //Создаем массивы главной
            //и побочной диагоналей
            double[] d = new double[n];
            double[] e = new double[n];

            //Заоплняем их значениями.
            Parallel.For(0, n - 1, i =>
            {
                d[i] = A[i, i];
                e[i] = A[i, i + 1];
            });
            d[n - 1] = A[n - 1, n - 1];

            //Вычисляем собственные числа и вектора трехдиагональной матрицы.
            bool result = smatrixevvTD(d, e, n, ref lambda, ref vector, epsilon, vectorsFound);

            return result;
        }
        #endregion

        #region smatrixevvAdd
        /// <summary>
        /// Метод вычисляет собственные числа и вектора матрицы, пополненной одним столбцом и строкой.
        /// Т.е. известны собственные числа и вектора исходной матрицы размерности N-1.
        /// </summary>
        /// <param name="A">Исходня матрица.</param>
        /// <param name="lambda">Собственные числа исходной матрицы (на входе - массив длины N-1, на выходе - массив длины N).</param>
        /// <param name="vector">Собственные числа исходной матрицы (на входе - матрица размерности N-1, на выходе - матрица размерности N).</param>
        /// <param name="epsilon">Точность вычисления.</param>
        /// <param name="vectorsFound">Индикатор, указывающий на надобность вычисления собственных векторов.
        /// <value>false</value> - если вычисление собственных векторов не требуется и <value>true</value> - в ином случае.</param>
        /// <returns></returns>
        public static bool smatrixevvAdd(double[,] A, ref double[] lambda, ref double[,] vector, double epsilon, bool vectorsFound)
        {
            //Получаем размерность матрицы
            int N = A.GetLength(0);

            //Создаем массивы по подобию трехдиагональной матрицы
            //главной и побочной диагоналей
            var D = new double[N];
            var E = new double[N - 1];
            //Выделяем последний (новый) столбец.
            var AN = new double[N];

            //Заполняем
            fixed(double* a = A, _AN = AN, d = D, l = lambda)
                vectorsAdd(_AN, a, l, d, N);            
            D[N - 1] = A[N - 1, N - 1];
            E[N - 2] = A[N - 2, N - 1];
            var tempAN = new double[N];

            //Вычисляем ограниченно-диагональную матрицу
            fixed (double* v = vector, _AN = AN, tA = tempAN)
                FormANkAdd(_AN, v, tA, N - 1);
            AN = tempAN;

            //Вычисляем собственные числа и вектора
            //ограниченно-диагональной матрицы.
            var L = new double[N];
            var Vt = new double[N, N];
            fixed (double* l = L, v = Vt, _d = D, _e = E, _a = AN)
                EigenValueAdd(_d, _e, _a, N, l, v, epsilon);

            //Вычисляем собственные вектора исходной
            //матрицы, если требуется.
            if (vectorsFound)
            {
                var tmp = new double[N, N];
                fixed(double* v = vector, t = tmp)
                    DiagonalVector(v, t, N);

                var V = (double[,])vector.Clone();
                vector = new double[N, N];
                fixed(double* v = V, _v = vector, t = tmp)
                    MatrixMultiply(t, v, _v, N);
            }

            lambda = L;

            return true;
        }
        #endregion          

    }

    public unsafe class smatrix2TD
    {
        const double machineepsilon = 5E-16;
        const double maxrealnumber = 1E300;
        const double minrealnumber = 1E-300;


        /*************************************************************************
       Преобразование симметричной матрицы, которая получена ее верхней или нижней
       треугольной частью к трехдиагональной матрице, используя простое ортогональное
       преобразование: Q'*A*Q=T.

       Входные параметры:
           A       -   Матрица для преобразования.
                       Массив с элементами [0..N-1, 0..N-1].
           N       -   Размерность матрицы A.
           IsUpper -   Формат хранения. Если IsUpper = True, то матрица А получена
                       ее верхней треугольной частью, нижняя треугольная часть не 
                       используется и не меняется алгоритмом, и наоборот. Если IsUpper 
                       = False.

       Выходные параметры:
           A       -   Матрицы T и Q в компактном виде (см. ниже).
           Tau     -   Массив факторов, которые формируют матрицу H(i).
                       Массив с элементами [0..N-2].
           D       -   Главная диагональ симметричной матрицы T.
                       Массив с элементами [0..N-1].
           E       -   Побочная диагональ симметричной матрицы T.
                       Массив с элементами [0..N-2].


         Если IsUpper=True, матрица Q представляет из себя продукт элементарных
         отражений

            Q = H(n-2) . . . H(2) H(0).

         Каждый H(i) имеет форму

            H(i) = I - tau * v * v'

         где tau - вещественный скаляр, и v вещественный вектор с
         v(i+1:n-1) = 0, v(i) = 1, v(0:i-1) хранится на выходе в
         A(0:i-1,i+1), и tau в TAU(i).

         Если IsUpper=False, матрица Q представляет из себя продукт элементарных
         отражений

            Q = H(0) H(2) . . . H(n-2).

         Каждый H(i) имеет форму

            H(i) = I - tau * v * v'

         где tau - вещественный скаляр, и v вещественный вектор с
         v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) зранящихся на выходе в A(i+2:n-1,i),
         и tau в TAU(i).

         Содержание A на выходе иллюстрируется следующим примером
         с n = 5:

         Если UPLO = 'U':                       Если UPLO = 'L':

           (  d   e   v1  v2  v3 )              (  d                  )
           (      d   e   v2  v3 )              (  e   d              )
           (          d   e   v3 )              (  v0  e   d          )
           (              d   e  )              (  v0  v1  e   d      )
           (                  d  )              (  v0  v1  v2  e   d  )

         где d и e обозначают диагональные и недиагональные элементы T, и vi
         обознаюает элементы вектора определения H(i).

       *************************************************************************/
        public static void smatrixtd(ref double[,] a, int n, bool isupper, ref double[] tau,
                                     ref double[] d, ref double[] e)
        {
            double alpha = 0.0d;
            double taui = 0.0d;
            double v = 0.0d;
            double[] t = new double[0];
            double[] t2 = new double[0];
            double[] t3 = new double[0];

            tau = new double[0];
            d = new double[0];
            e = new double[0];

            if (n <= 0)
            {
                return;
            }

            t = new double[n + 1];
            t2 = new double[n + 1];
            t3 = new double[n + 1];
            d = new double[n];

            if (n > 1)
            {
                tau = new double[n - 1];
            }

            if (n > 1)
            {
                e = new double[n - 1];
            }

            fixed (double* _a = a, _t = t, _e = e, _tau = tau, _d = d)
                if (isupper)
                {

                    //
                    // Уменьшаем верхний треугольник A
                    //
                    int i = n - 1;
                    for (double* stj = _a + n - 1, spj = _a + 1, spi = _a + (n - 2) * n - 1, E = _e + n - 2, D = _d + n - 1, Tau = _tau + n - 2;
                    stj >= spj;
                    stj--, spi -= (n + 1), i--, E--, D--, Tau--)
                    {

                        //
                        // Генерируем элементарное отражение H() = E - tau * v * v'
                        //
                        if (stj >= _a + 2)
                        {
                            for (double* sti = stj, ti = _t + 2; sti <= spi; sti += n, ti++)
                            {
                                *ti = *sti;
                            }
                        }
                        _t[1] = *(spi + n);
                        generatereflection(_t, i, out taui);

                        if (stj >= _a + 2)
                        {
                            //for (i_ = 0; i_ <= i - 1; i_++)
                            for (double* sti = stj, ti = _t + 2; sti <= spi; sti += n, ti++)
                            {
                                *sti = *ti;
                                //a[i_, i + 1] = t[i_ + 2];
                            }
                        }
                        *(spi + n) = _t[1];
                        *E = *(spi + n);
                        fixed (double* _t3 = t3, _t2 = t2)
                            if (taui != 0.0d)
                            {

                                //
                                // Применяем отражение H с двух сторон к A
                                //
                                *(spi + n) = 1.0d;

                                //
                                // Рассчитываем  x := tau * A * v  storing x in TAU
                                //
                                for (double* sti = stj, ti = _t + 1; sti <= spi + n; sti += n, ti++)
                                {
                                    *ti = *sti;
                                }
                                symmetricmatrixvectormultiply(_a, isupper, 0, i - 1, t, taui, _t3, n);
                                for (double* st = _tau, sp = _tau + i - 1, t3i = _t3 + 1; st <= sp; st++, t3i++)
                                {
                                    *st = *t3i;
                                }

                                //
                                // Рассчитываем  w := x - 1/2 * tau * (x'*v) * v
                                //
                                v = 0.0d;
                                for (double* sti = stj, _taui = _tau; sti <= spi + n; sti += n, _taui++)
                                {
                                    v += *_taui * *sti;
                                }
                                alpha = -(0.5d * taui * v);
                                for (double* sti = stj, _taui = _tau; sti <= spi + n; sti += n, _taui++)
                                {
                                    *_taui += alpha * *sti;
                                }

                                //
                                // Применяем преобразование как обновление rank-2:
                                //    A := A - v * w' - w * v'
                                //
                                for (double* sti = stj, _ti = _t + 1; sti <= spi + n; sti += n, _ti++)
                                {
                                    *_ti = *sti;
                                }
                                for (double* st = _tau, _t3i = _t3 + 1, sp = _tau + i; st <= sp; st++, _t3i++)
                                {
                                    *_t3i = *st;
                                }
                                symmetricrank2update(_a, isupper, 0, i - 1, t, t3, _t2, -1.0d, n);
                                *(spi + n) = *E;
                            }
                        *D = *(spi + n + n);
                        *Tau = taui;
                    }
                    *_d = *_a;
                }
                else
                {

                    //
                    // Уменьшаем нижний треугольник A
                    //
                    int i = 0;
                    for (double* stj = _a + n, spj = _a + n * n - 2, spi = _a + (n - 1) * n, E = _e, D = _d, Tau = _tau;
                    stj <= spj;
                    stj += (n + 1), spi++, i++, E++, D++, Tau++)
                    {

                        //
                        // Генерируем элементарное отражение H = E - tau * v * v'
                        //
                        for (double* sti = stj, ti = _t + 1; sti <= spi; sti += n, ti++)
                        {
                            *ti = *sti;
                            //t[i_] = a[i_ + i, i];
                        }
                        generatereflection(_t, n - i - 1, out taui);
                        for (double* sti = stj, ti = _t + 1; sti <= spi; sti += n, ti++)
                        {
                            *sti = *ti;
                        }
                        *E = *stj;
                        fixed (double* _t3 = t3, _t2 = t2)
                            if (taui != 0.0d)
                            {

                                //
                                // Применяем отражение H с двух сторон к A
                                //
                                *stj = 1.0d;

                                //
                                // Рассчитываем  x := tau * A * v  storing y in TAU
                                //
                                for (double* sti = stj, ti = _t + 1; sti <= spi; sti += n, ti++)
                                {
                                    *ti = *sti;
                                }
                                symmetricmatrixvectormultiply(_a, isupper, i + 1, n - 1, t, taui, _t2, n);
                                for (double* st = Tau + i, sp = Tau + n - 2, t2i = _t2; st <= sp; st++, t2i++)
                                {
                                    *st = *t2i;
                                }

                                //
                                // Рассчитываем w := x - 1/2 * tau * (x'*v) * v
                                //
                                v = 0.0d;
                                for (double* st = stj, _taui = Tau + i; st <= spi; st += n, _taui++)
                                {
                                    v += *_taui * *st;
                                }
                                alpha = -(0.5 * taui * v);
                                for (double* st = stj, _taui = Tau + i; st <= spi; st += n, _taui++)
                                {
                                    *_taui += alpha * *st;
                                }

                                //
                                // Применяем преобразование как обновление rank-2:
                                //     A := A - v * w' - w * v'
                                //
                                // 
                                for (double* sti = stj, ti = _t + 1; sti <= spi; sti += n, ti++)
                                {
                                    *ti = *sti;
                                }
                                for (double* st = Tau + i, t2i = _t2 + 1, sp = Tau + n - 2; st <= sp; st++, t2i++)
                                {
                                    *t2i = *st;
                                }
                                symmetricrank2update(_a, isupper, i + 1, n - 1, t, t2, _t3, -1.0d, n);
                                *stj = *E;
                            }
                        *D = *(stj - n);
                        *Tau = taui;
                    }
                    _d[n - 1] = _a[n * n - 1];
                }
        }

        /*************************************************************************
      Генерирует элементарное преобразование отражения

      Подпрограмма генерирует элементарное отражение H порядка N, так, что для
      данного X, имеет место равенство:

          ( X(1) )   ( Beta )
      H * (  ..  ) = (  0   )
          ( X(n) )   (  0   )

      где
                    ( V(1) )
      H = 1 - Tau * (  ..  ) * ( V(1), ..., V(n) )
                    ( V(n) )

      где первая компонента вектора V равна 1.

      Входные параметры:
          X   -   Вектор. Массив с индексами в диапазоне [1..N].
          N   -   Порядок отражения.

      ВЫходные параметры:
          X   -   Кмопоненты от 2 до N, размещенные в векторе V.
                  Первая компонента заменяется с параметром Beta.
          Tau -   Скалярное значение Tau. Если X нулевой вектор, Tau равно 0,
                  иначе 1 <= Tau <= 2.

      Эта подпрограмма является модификацией DLARFG подпрограммы из
      библиотеки LAPACK.

      *************************************************************************/
        static void generatereflection(double* x, int n, out double tau)
        {
            double alpha = 0.0d;
            double xnorm = 0.0d;
            double v = 0.0d;
            double beta = 0.0d;
            double mx = 0.0d;
            double s = 0.0d;

            tau = 0.0d;

            if (n <= 1)
                return;

            //
            // Если требуется, масштабируем (во избежание переполнения/недополнения
            // во время промежуточных вычислений).
            //
            for (double* st = x + 1, sp = x + n; st <= sp; st++)
            {
                mx = Math.Max(Math.Abs(*st), mx);
            }
            s = 1.0d;
            if (mx != 0.0d)
            {
                if (mx <= (minrealnumber / machineepsilon))
                {
                    s = minrealnumber / machineepsilon;
                    v = 1 / s;
                    for (double* st = x + 1, sp = x + n; st <= sp; st++)
                    {
                        *st = v * *st;
                    }
                    mx = mx * v;
                }
                else
                {
                    if (mx >= (maxrealnumber * machineepsilon))
                    {
                        s = maxrealnumber * machineepsilon;
                        v = 1.0d / s;
                        for (double* st = x + 1, sp = x + n; st <= sp; st++)
                        {
                            *st = v * *st;
                        }
                        mx = mx * v;
                    }
                }
            }

            //
            // XNORM = DNRM2( N-1, X, INCX )
            //
            alpha = x[1];
            if (mx != 0.0d)
            {
                for (double* st = x + 2, sp = x + n; st <= sp; st++)
                {
                    xnorm += Sqr(*st / mx);
                }
                xnorm = Math.Sqrt(xnorm) * mx;
            }
            if (xnorm == 0.0d)
            {

                //
                // H  =  I
                //
                x[1] = x[1] * s;
                return;
            }

            //
            // Основной случай
            //
            mx = Math.Max(Math.Abs(alpha), Math.Abs(xnorm));
            beta = -(mx * Math.Sqrt(Sqr(alpha / mx) + Sqr(xnorm / mx)));
            if (alpha < 0.0d)
            {
                beta = -beta;
            }
            tau = (beta - alpha) / beta;
            v = 1.0d / (alpha - beta);
            for (double* st = x + 2, sp = x + n; st <= sp; st++)
            {
                *st = v * *st;
            }
            x[1] = beta;

            //
            // Обратное масштабирование выхода
            //
            x[1] = x[1] * s;
        }

        static void symmetricmatrixvectormultiply(double* a, bool isupper, int i1, int i2,
                                                         double[] X, double alpha, double* y, int N)
        {
            int n = i2 - i1 + 1;

            if (n <= 0)
            {
                return;
            }

            //
            // Примем A = L + D + U, где
            //  L - строго нижний треугольник (нулевая главная диагональ)
            //  D - диагональ
            //  U - строго верхний треугольник (нулевая главная диагональ)
            //
            // A*x = L*x + D*x + U*x
            //
            // Calculate D*x first
            //          
            fixed (double* x = X)
                for (double* st = a + N * i1 + i1, sp = a + N * i2 + i2, _x = x + 1, _y = y + 1; st <= sp; st += N + 1, _x++, _y++)
                {
                    *_y = *st * *_x;
                }

            //
            // Добавляем L*x + U*x
            //
            fixed (double* x = X)
                if (isupper)
                {
                    for (int i = i1; i < i2; i++)
                    {

                        //
                        // Добавляем L*x к результату
                        //
                        double v = x[i - i1 + 1];
                        for (double* st = a + i * N + i + 1, sp = st + n - i - 2 + i1, _y = y + i - i1 + 2; st <= sp; st++, _y++)
                        {
                            *_y = *_y + v * *st;
                        }

                        //
                        // Добавляем U*x к результату
                        //
                        v = 0.0d;
                        for (double* st = a + i * N + i + 1, sp = st + n - i - 2 + i1, _x = x + i - i1 + 2; st <= sp; st++, _x++)
                        {
                            v += *_x * *st;
                        }
                        y[i - i1 + 1] += v;
                    }
                }
                else
                {
                    for (int i = i1 + 1; i <= i2; i++)
                    {

                        //
                        // Добавляем L*x к результату
                        //
                        double v = 0.0d;
                        for (double* _x = x + 1, st = a + i * N + i1, sp = st + i - 1 - i1; st <= sp; st++, _x++)
                        {
                            v += *_x * *st;
                        }
                        y[i - i1 + 1] = y[i - i1 + 1] + v;

                        //
                        // Добавляем U*x к результату
                        //
                        v = x[i - i1 + 1];
                        for (double* _y = y + 1, st = a + i * N + i1, sp = st + i - 1 - i1; st <= sp; st++, _y++)
                        {
                            *_y = *_y + v * *st;
                        }
                    }
                }

            for (double* st = y + 1, sp = y + n; st <= sp; st++)
            {
                *st = alpha * *st;
            }
        }

        static void symmetricrank2update(double* a, bool isupper, int i1, int i2,
                                                double[] X, double[] Y, double* t, double alpha, int N)
        {
            int i = 0;

            fixed (double* x = X, y = Y)
                if (isupper)
                {
                    for (i = i1; i <= i2; i++)
                    {
                        double v = x[i + 1 - i1];
                        for (double* st = y + i + 1 - i1, _t = t + i + 1 - i1, sp = y + i2 - i1 + 1; st <= sp; st++, _t++)
                        {
                            *_t = v * *st;
                        }

                        v = y[i + 1 - i1];
                        for (double* st = x + i + 1 - i1, _t = t + i + 1 - i1, sp = x + i2 - i1 + 1; st <= sp; st++, _t++)
                        {
                            *_t = *_t + v * *st;
                        }

                        for (double* st = t + i + 1 - i1, sp = t + i2 - i1 + 1; st <= sp; st++)
                        {
                            *st = alpha * *st;
                        }

                        for (double* st = a + N * i + i, sp = st + i2 - i, _t = t + i + 1 - i1; st <= sp; _t++, st++)
                        {
                            *st += *_t;
                        }
                    }
                }
                else
                {
                    for (i = i1; i <= i2; i++)
                    {
                        double v = x[i + 1 - i1];

                        for (double* st = y + 1, _t = t + 1, sp = y + i - i1 + 1; st <= sp; st++, _t++)
                        {
                            *_t = v * *st;
                        }

                        v = y[i + 1 - i1];
                        for (double* st = x + 1, _t = t + 1, sp = x + i - i1 + 1; st <= sp; st++, _t++)
                        {
                            *_t += v * *st;
                        }

                        for (double* st = t + 1, sp = t + i - i1 + 1; st <= sp; st++)
                        {
                            *st = alpha * *st;
                        }

                        for (double* st = a + N * i + i1, sp = a + N * i + i, _t = t + 1; st <= sp; st++, _t++)
                        {
                            *st += *_t;
                        }
                    }
                }
        }



        /*************************************************************************
       Распаковка симметричной матрицы Q из трехдиагональной формы.

       Входные параметры:
           A       -   Результат подпрограммы SMatrixTD.
           N       -   Размерность матрицы A.
           IsUpper -   Формат хранения (параметр подпрограммы SMatrixTD).
           Tau     -   результат подпрограммы SMatrixTD.

       Выходные параметры:
           Q       -   Матрица преобразования.
                       Массив с элементами. [0..N-1, 0..N-1].

       *************************************************************************/
        public static void smatrixtdunpackq(double[,] a, int n, bool isupper, double[] tau, ref double[,] q)
        {
            double[] v = new double[0];
            double[] work = new double[0];

            q = new double[n, n];

            if (n == 0)
            {
                return;
            }

            //
            // Инициализация
            //

            v = new double[n + 1];
            work = new double[n];

            fixed (double* _q = q)
                for (double* st = _q, sp = _q + n * n - 1; st <= sp; st += n + 1)
                {
                    *st = 1.0d;
                }

            //
            // Распаковка Q
            //
            fixed (double* w = work, _q = q, _v = v, _a = a)
                if (isupper)
                {
                    for (int j = 0; j < n - 1; j++)
                    {

                        //
                        // Применяем H(i)
                        //
                        for (int i = 0; i < j; i++)
                        {
                            _v[i + 1] = _a[i * n + j + 1];
                        }
                        _v[j + 1] = 1.0d;
                        applyreflectionfromtheleft(_q, tau[j], _v, 0, j, 0, n - 1, w);
                    }
                }
                else
                {
                    for (int j = n - 2; j >= 0; j--)
                    {

                        //
                        // Применяем H(i)
                        //
                        for (int i = 1; i <= n - j - 1; i++)
                        {
                            _v[i] = _a[(i + j) * n + j];
                        }
                        _v[1] = 1.0d;
                        applyreflectionfromtheleft(_q, tau[j], _v, j + 1, n - 1, 0, n - 1, w);
                    }
                }
        }

        /*************************************************************************
        Применение элементарного отражения к прямоугольной матрице размераe MxN.

        Алгоритм предварительно умножает матрицу на преобразование элементарного
        отражения, который получен столбцом V и скаляром Tau (см. описание процедуры
        GenerateReflection). Преобразуется лишь часть матрицы (строки от M1 до M2,
        столбцы от N1 до N2), а не вся матрица. Меняются только элементы данной
        подматрицы.

        Входные параметры:
            C       -   Матрица для преобразования.
            Tau     -   Скаляр определяющий преобразование.
            V       -   Столбец, определяющий преобразование.
                        Массив с индексами в диапазоне [1..M2-M1+1].
            M1, M2  -   Диапазон строк для преобразования.
            N1, N2  -   Диапазон столбцов для преобразования.
            WORK    -   Рабочий массив, индексы которого лежат в диапазоне от N1 до N2.

        Выходные параметры:
            C       -   Результат умножения входной матрицы C с помощью преобразующей
                        матрицы полученной из Tau и V.
                        Если N1>N2 или M1>M2, C не меняется.

        *************************************************************************/
        static void applyreflectionfromtheleft(double* c, double tau, double* v, int m1,
                                                      int m2, int n1, int n2, double* work)
        {
            double t = 0.0d;

            if ((tau == 0.0d || n1 > n2) || m1 > m2)
            {
                return;
            }

            //
            // w := C' * v
            //
            for (double* st = work + n1, sp = work + n2; st <= sp; st++)
            {
                *st = 0.0d;
            }


            for (double* sti = c + (n2 + 1) * m1, spi = c + (n2 + 1) * m2, _v = v + 1; sti <= spi; sti += n2 + 1, _v++)
            {
                for (double* stj = sti + n1, spj = sti + n2, w = work + n1; stj <= spj; stj++, w++)
                {
                    *w += *_v * *stj;
                }
            }

            //
            // C := C - tau * v * w'
            //
            for (double* _v = v + 1, sti = c + (n2 + 1) * m1, spi = c + (n2 + 1) * m2; sti <= spi; sti += n2 + 1, _v++)
            {
                t = *_v * tau;
                for (double* stj = sti + n1, spj = sti + n2, w = work + n1; stj <= spj; stj++, w++)
                {
                    *stj = *stj - t * *w;
                }
            }
        }
        
        static double Sqr(double p)
        {
            return p * p;
        }

    }  
    

    [StructLayout(LayoutKind.Explicit)]
    public struct BigDouble
    {
        //Число от 1.0 до 9.99
        [FieldOffset(0)]
        double x;
        //Степень числа
        [FieldOffset(8)]
        int y;

        
        public double X
        {
            get { return x; }
            set { x = value; }
        }        
        public int Y
        {
            get { return y; }
            set { y = value; }
        }

        /// <summary>
        /// Конструктор, создающий число типа BigDouble вида x*10^y 
        /// из строкового представления типа Double вида XEY.
        /// </summary>
        /// <param name="S"></param>
        public BigDouble(string S)
        {
            string[] sa = S.Split(new char[] { 'E' });

            this.x = Double.Parse(sa[0]);
            if (sa.Length > 1)
                this.y = int.Parse(sa[1]);
            else
                this.y = 0;

            while (Math.Abs(this.X) > 9.99d)
            {
                this.X /= 10.0d;
                this.Y += 1;
            }

            if (this.X != 0.0d)
                while (Math.Abs(this.X) < 1.0d)
                {
                    this.X *= 10.0d;
                    this.Y -= 1;
                }
        }

        /// <summary>
        /// Конструктор, создающий число типа BigDouble 
        /// вида x*10^y по заданным значениям.
        /// </summary>
        /// <param name="d">Число от 1.0 до 9.99</param>
        /// <param name="n">Степень числа</param>
        public BigDouble(double d, int n)
        {
            this.x = d;
            this.y = n;

            while (Math.Abs(this.X) > 9.99d)
            {
                this.X /= 10.0d;
                this.Y += 1;
            }

            if (this.X != 0.0d)
                while (Math.Abs(this.X) < 1.0d)
                {
                    this.X *= 10.0d;
                    this.Y -= 1;
                }
        }


        //Перегрузка соотвествующих операторов

        public static BigDouble operator +(BigDouble d1, BigDouble d2)
        {
            if (d1.Y > d2.Y)
            {
                for (int i = 0; i < d1.Y - d2.Y; i++)
                {
                    d1.X *= 10.0d;
                }

                d1.Y = d2.Y;
            }
            else if (d2.Y > d1.Y)
            {
                for (int i = 0; i < d2.Y - d1.Y; i++)
                {
                    d2.X *= 10.0d;
                }

                d2.Y = d1.Y;
            }

            return Sum(d1, d2);
        }

        public static BigDouble operator -(BigDouble d1, BigDouble d2)
        {
            if (d1.Y > d2.Y)
            {
                for (int i = 0; i < d1.Y - d2.Y; i++)
                {
                    d1.X *= 10.0d;
                }

                d1.Y = d2.Y;
            }
            else if (d2.Y > d1.Y)
            {
                for (int i = 0; i < d2.Y - d1.Y; i++)
                {
                    d2.X *= 10.0d;
                }

                d2.Y = d1.Y;
            }

            return Minus(d1, d2);
        }

        public static BigDouble operator *(BigDouble d1, BigDouble d2)
        {
            d1.X *= d2.X;
            d1.Y += d2.Y;

            while(Math.Abs(d1.X) > 9.99d)
            {
                d1.X /= 10.0d;
                d1.Y += 1;
            }

            if (d1.X != 0.0d)
                while (Math.Abs(d1.X) < 1.0d)
                {
                    d1.X *= 10.0d;
                    d1.Y -= 1;
                }

            return d1;
        }

        public static BigDouble operator /(BigDouble d1, BigDouble d2)
        {
            d1.X /= d2.X;
            d1.Y -= d2.Y;

            while(Math.Abs(d1.X) > 9.99d)
            {
                d1.X /= 10.0d;
                d1.Y += 1;
            }

            if (d1.X != 0.0d)
                while (Math.Abs(d1.X) < 1.0d)
                {
                    d1.X *= 10.0d;
                    d1.Y -= 1;
                }

            return d1;
        }

        public static explicit operator BigDouble(double d)
        {
            return new BigDouble(d.ToString());
        }

        public static explicit operator double(BigDouble d)
        {
            if (d.Y > 307)
                return Double.PositiveInfinity;
            else if (d.Y < -323)
                return 0.0d;

            return d.X * Math.Pow(10.0d, (double)d.Y);
        }

        public static explicit operator string(BigDouble d)
        {
            return d.ToString();
        }


        /// <summary>
        /// Метод вычисляет сумму двух чисел типа BigDouble одной степени.
        /// </summary>
        /// <param name="d1">Первое число.</param>
        /// <param name="d2">Второе число.</param>
        /// <returns>Сумма двух чисел.</returns>
        static BigDouble Sum(BigDouble d1, BigDouble d2)
        {  
            d1.X += d2.X;
            while(Math.Abs(d1.X) > 9.99d)
            {
                d1.X /= 10.0d;
                d1.Y += 1;
            }

            if (d1.X != 0.0d)
                while (Math.Abs(d1.X) < 1.0d)
                {
                    d1.X *= 10.0d;
                    d1.Y -= 1;
                }

            return d1;
        }

        /// <summary>
        /// Метод вычисляет разность двух чисел типа BigDouble одной степени.
        /// </summary>
        /// <param name="d1">Первое число.</param>
        /// <param name="d2">Второе число.</param>
        /// <returns>Разность двух чисел.</returns>
        static BigDouble Minus(BigDouble d1, BigDouble d2)
        {
            d1.X -= d2.X;
            while(Math.Abs(d1.X) > 9.99d)
            {
                d1.X /= 10.0d;
                d1.Y += 1;
            }

            if (d1.X != 0.0d)
                while (Math.Abs(d1.X) < 1.0d)
                {
                    d1.X *= 10.0d;
                    d1.Y -= 1;
                }

            return d1;           
        }

        /// <summary>
        /// Метод возводит заданное число типа BigDouble в указанную степень.
        /// </summary>
        /// <param name="d">Число.</param>
        /// <param name="n">Степень возведения.</param>
        /// <returns>Число возведенное в заданную степень.</returns>
        public static BigDouble Pow(BigDouble d, int n)
        {
            d.Y *= n;
            if (n == 2)
                d.X *= d.X;
            else
            {
                for (int i = 0; i < n; i++)
                    d.X *= d.X;
            }

            while (Math.Abs(d.X) > 9.99d)
            {
                d.X /= 10.0d;
                d.Y += 1;
            }

            if (d.X != 0.0d)
                while (Math.Abs(d.X) < 1.0d)
                {
                    d.X *= 10.0d;
                    d.Y -= 1;
                }

            return d;
        }

        /// <summary>
        /// Метод вычисляет квадратный корень числа типа BigDouble.
        /// </summary>
        /// <param name="d">Заданное число.</param>
        /// <returns>Корень квадратный из заданного числа.</returns>
        public static BigDouble Sqrt(BigDouble d)
        {
            if (d.Y % 2 == 0)
                d.Y /= 2;
            else
            {
                d.Y--;
                d.Y /= 2;
                d.X *= 10.0d;
            }

            d.X = Math.Sqrt(d.X);

            while (Math.Abs(d.X) > 9.99d)
            {
                d.X /= 10.0d;
                d.Y += 1;
            }
            
            if (d.X != 0.0d)
                while (Math.Abs(d.X) < 1.0d)
                {
                    d.X *= 10.0d;
                    d.Y -= 1;
                }

            return d;
        }

        /// <summary>
        /// Метод возвращает строковое представление типа BigDouble.
        /// </summary>
        /// <returns>Строковое представление числа.</returns>
        public override string ToString()
        {
            if (this.Y != 0)
                return String.Format("{0}E{1}", this.X, this.Y);

            return this.X.ToString();
        }
    }
}

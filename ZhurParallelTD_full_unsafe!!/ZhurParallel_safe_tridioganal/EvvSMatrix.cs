using System;
using System.Threading.Tasks;

namespace ZhurParallelTDusF
{
    public unsafe partial class EigenValueAndVectorProblem
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

            if (L2 - L1 == 0.0d)
                return L1;

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

            while (N < 0)
            {
                epsilon *= epsilon;

                _N = Math.Log((L2 - L1) / epsilon);
                N = (int)Math.Ceiling(_N * 3.32192809488736d);
            }
            
            
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


        #region smatrixevv
        /// <summary>
        /// Метод вычисляет собственные вектора и собственные значения исходной симметричной матрицы A 
        /// c заданной точностью epsilon, путем преобразования ее к трехдиагональному виду.
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
    }
}

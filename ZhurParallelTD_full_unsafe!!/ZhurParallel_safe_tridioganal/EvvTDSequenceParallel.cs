using System;
using System.Threading.Tasks;

namespace ZhurParallelTDusF
{
    public unsafe partial class EigenValueAndVectorProblem
    {
        #region EigenValueAndVector2x2TD
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
        #endregion

        #region EigenValue
        /// <summary>
        /// Метод рассчитывает собственные вектора и числа заданной ограниченнно-диагональной матрицы с указанной точностью.
        /// </summary>
        /// <param name="D">Указатель на массив, составляющий главную диагональ исходной матрицы.</param>
        /// <param name="A">Указатель на массив, составляющий ненулевой столбец ограниченно-диагональной матрицы.</param>
        /// <param name="n">Размерность ограниченно-диагональной матрицы (или шаг вычисления).</param>
        /// <param name="lambda">Указатель на массив, составляющий сосбтвенные числа.</param>
        /// <param name="vector">Указатель на массив, составляющий последнюю строку матрицы собственных векторов.</param>
        /// <param name="epsilon">Точность вычисления.</param>
        static void EigenValue(double* D, double* A, int n, double* lambda, double* vector, double epsilon, bool up)
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
            EigenVectorsLD(vector, lambda, A, D, n, up);
        }
        #endregion

        #region EigenVectorsLD
        /// <summary>
        /// Метод вычисляет собственные вектора исходной ограниченно-диагональной матрицы.
        /// </summary>
        /// <param name="V">Указатель на последнюю строку матрицы собственных векторов ограниченно-диагональной матрицы.</param>
        /// <param name="lambda">Указатель на массив собственных чисел ограниченно-диагональной матрицы.</param>
        /// <param name="A">Указатель на ненулевой столбец ограниченно-диагональной матрицы.</param>
        /// <param name="D">Указатель на массив, составляющий главную диагональ матрицы.</param>
        /// <param name="n">Размерность матрицы</param>
        static void EigenVectorsLD(double* V, double* lambda, double* A, double* D, int n, bool up)
        {
            //Вычисляем собственные вектора и
            //заносим последнюю или первую строку в массив.
            Parallel.For(0, n, j =>
            {
                double temp = 1.0d;

                for (double* i = D, _a = A, i_end = D + n - 1; i < i_end; i++, _a++)
                {

                    var tt = (lambda[j] - *i);
                    double t;
                    if (tt == 0)
                        t = 0;
                    else
                        t = Sqr(*_a / tt);

                    temp += t;
                }

                if (up)
                {
                    var del = lambda[j] - *D;
                    if (del == 0)
                        V[j] = 0;
                    else
                        V[j] = *A / (del * Math.Sqrt(temp));
                }
                else
                    V[j] = Math.Sqrt(1.0d / temp);
            });

        }
        #endregion

        #region EigenVectors
        /// <summary>
        /// Находит собственные вектора симметричной трехдиагональной матрицы по известным собственным значениям.
        /// Используется метод прогонки.
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
                fixed (double* l = L, v = V, _a = a, _d = d)
                    EigenValue(_d, _a, i + 1, l, v, epsilon, false);
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
        /// <param name="A">Трехдиагональная матрица, представленная двумерным массивом.</param>        
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

    }
}

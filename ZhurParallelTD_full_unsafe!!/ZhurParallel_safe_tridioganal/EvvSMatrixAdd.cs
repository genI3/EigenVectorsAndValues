using System;
using System.Threading.Tasks;

namespace ZhurParallelTDusF
{
    public unsafe partial class EigenValueAndVectorProblem
    {
        #region FormANkAdd
        /// <summary>
        /// Заполняем массивы главной диагонали, составляющей собственные числа, и
        /// последнего (нового) столбца матрицы.
        /// </summary>
        /// <param name="AN">Указатель на массив, составляющий последний (новый) столбец матрицы.</param>
        /// <param name="A">Указатель на исходную матрицу.</param>
        /// <param name="lambda">Указатель на массив собственных чисел.</param>
        /// <param name="D">Указатель на массив, составляющий главную диагональ.</param>
        /// <param name="n">Размерность матрицы.</param>
        static void FormANkAdd(double* AN, double* A, double* lambda, double* D, int n, double* v)
        {
            //Заполняем массив последнего столбца
            Parallel.For(0, n - 1, i =>
            {
                AN[i] = A[i * n + n - 1];
            });

            //Заполняем главную диагональ.
            Parallel.For(0, n - 1, i =>
            {
                D[i] = lambda[i];
            });

            //Вычисляем ограниченно-диагональную матрицу
            FormAColumnAdd(AN, v, n - 1);
        }

        /// <summary>
        /// Формирует ограниченно-диагональную матрицу ANk для последнего шага.
        /// </summary>
        /// <param name="AN">Указатель на массив, составляющий последний (новый) столбец исходной матрицы.</param>
        /// <param name="V">Указатель на матрицу собственных векторов матрицы размерности N-1, вырезанной из исходной.</param>
        /// <param name="tempA">Указатель на временный массив.</param>
        /// <param name="n">Размерность матрицы.</param>
        static void FormAColumnAdd(double* AN, double* V, int n)
        {
            var tempA = new double[n - 1];

            fixed (double* ta = tempA)
            {
                var tA = ta;

                Parallel.For(0, n, i =>
                {
                    var temp = 0.0d;
                    for (int m = 0; m < n; m++)
                        temp += AN[m] * V[m * n + i];
                    tA[i] = temp;
                });

                Parallel.For(0, n, i =>
                {
                    AN[i] = tA[i];
                });
            }
        }
        #endregion

        #region EigenValueAdd
        /// <summary>
        /// Метод рассчитывает собственные вектора и числа заданной ограниченнно-диагональной матрицы с указанной точностью.
        /// </summary>
        /// <param name="D">Указатель на массив, составляющий главную диагональ исходной матрицы.</param>
        /// <param name="A">Указатель на массив, составляющий ненулевой столбец ограниченно-диагональной матрицы.</param>
        /// <param name="n">Размерность ограниченно-диагональной матрицы (или шаг вычисления).</param>
        /// <param name="n">Размерность исходной матрицы</param>
        /// <param name="lambda">Указатель на массив, составляющий сосбтвенные числа.</param>
        /// <param name="vector">Указатель на массив, составляющий последнюю строку матрицы собственных векторов.</param>
        /// <param name="epsilon">Точность вычисления.</param>
        static void EigenValueAdd(double* D, double* A, int n, double* lambda, double* vector, double epsilon)
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
            EigenVectorAdd(vector, lambda, A, D, n);
        }
        #endregion

        #region EigenVectorAdd
        /// <summary>
        /// Метод вычисляет собственные вектора исходной ограниченно-диагональной матрицы.
        /// </summary>
        /// <param name="V">Указатель на матрицу собственных векторов ограниченно-диагональной матрицы.</param>
        /// <param name="lambda">Указатель на массив собственных чисел ограниченно-диагональной матрицы.</param>
        /// <param name="A">Указатель на ненулевой столбец ограниченно-диагональной матрицы.</param>
        /// <param name="D">Указатель на массив, составляющий главную диагональ матрицы.</param>
        /// <param name="n">Размерность матрицы</param>
        static void EigenVectorAdd(double* V, double* lambda, double* A, double* D, int n)
        {
            //Вычисляем собственные вектора.
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

            //Создаем массив по подобию трехдиагональной 
            //матрицы главной диагонали
            var D = new double[N];
            //Выделяем последний (новый) столбец.
            var AN = new double[N - 1];

            //Заполняем
            fixed (double* a = A, _AN = AN, d = D, l = lambda, v = vector)
                FormANkAdd(_AN, a, l, d, N, v);
            D[N - 1] = A[N - 1, N - 1];

            //Вычисляем собственные числа и вектора
            //ограниченно-диагональной матрицы.
            var L = new double[N];
            var Vt = new double[N, N];
            fixed (double* l = L, v = Vt, _d = D, _a = AN)
                EigenValueAdd(_d, _a, N, l, v, epsilon);

            //Вычисляем собственные вектора исходной
            //матрицы, если требуется.
            if (vectorsFound)
            {
                var tmp = new double[N, N];
                fixed (double* v = vector, t = tmp)
                    DiagonalVector(v, t, N);

                vector = new double[N, N];
                fixed (double* v = Vt, _v = vector, t = tmp)
                    MatrixMultiply(t, v, _v, N);
            }

            lambda = L;

            return true;
        }
        #endregion
    }
}

using System;
using System.Threading.Tasks;

namespace ZhurParallelTDusF
{
    public unsafe partial class EigenValueAndVectorProblem
    {
        #region EigenValueAndVectorXxY
        /// <summary>
        /// 
        /// </summary>
        /// <param name="d"></param>
        /// <param name="e"></param>
        /// <param name="n1"></param>
        /// <param name="n2"></param>
        /// <param name="oldL1"></param>
        /// <param name="oldV1"></param>
        /// <param name="oldL2"></param>
        /// <param name="oldV2"></param>
        /// <param name="newL"></param>
        /// <param name="newV"></param>
        /// <param name="epsilon"></param>
        static void EigenValueAndVectorXxY(double d, double e1, double e2, int n1, int n2, double* oldL1, double* oldV1,
                                           double* oldL2, double* oldV2, double* newL, double epsilon)
        {
            var AN = new double[n1 + n2];
            var D = new double[n1 + n2 + 1];

            var tempL = new double[n1 + n2];

            fixed (double* l = tempL)
                FormL(l, oldL1, oldL2, n1, n2);

            fixed (double* a = AN)
                FormANkXxY(e1, e2, a, oldV1, oldV2, n1, n2);

            Array.Sort(tempL, AN);

            fixed (double* l = tempL, _d = D)
                FormDkXxY(_d, l, n1 + n2 + 1, d);

            fixed (double* _d = D, _a = AN)
                EigenValueNew(_d, _a, n1 + n2 + 1, newL, epsilon);

        }
        #endregion

        #region EigenValueAndVectorXxX
        /// <summary>
        /// 
        /// </summary>
        /// <param name="d"></param>
        /// <param name="e"></param>
        /// <param name="start_poz"></param>
        /// <param name="n"></param>
        /// <param name="oldL1"></param>
        /// <param name="oldV1"></param>
        /// <param name="oldL2"></param>
        /// <param name="oldV2"></param>
        /// <param name="newL"></param>
        /// <param name="newV"></param>
        /// <param name="epsilon"></param>
        /// <param name="up"></param>
        static void EigenValueAndVectorXxX(double d, double e1, double e2, int start_poz, int n, double* oldL1, double* oldV1,
                                           double* oldL2, double* oldV2, double* newL, double epsilon)
        {
            int n1 = (n - 1) / 2;

            var AN = new double[n - 1];
            var D = new double[n];

            var tempL = new double[n - 1];

            fixed (double* l = tempL)
                FormL(l, oldL1, oldL2, n1, n1);

            fixed (double* a = AN)
                FormANkXxY(e1, e2, a, oldV1, oldV2, n1, n1);

            Array.Sort(tempL, AN);

            fixed (double* l = tempL, _d = D)
                FormDkXxY(_d, l, n, d);

            fixed (double* _d = D, _a = AN)
                EigenValueNew(_d, _a, n, newL, epsilon);
        }
        #endregion

        #region FormL
        /// <summary>
        /// 
        /// </summary>
        /// <param name="L"></param>
        /// <param name="L1"></param>
        /// <param name="L2"></param>
        /// <param name="n1"></param>
        /// <param name="n2"></param>
        static void FormL(double* L, double* L1, double* L2, int n1, int n2)
        {
            Parallel.For(0, n1, i =>
            {
                L[i] = L1[i];
            });

            Parallel.For(n1, n1 + n2, i =>
            {
                L[i] = L2[i - n1];
            });
        }
        #endregion

        #region FormANkXxY
        /// <summary>
        /// 
        /// </summary>
        /// <param name="e1"></param>
        /// <param name="e2"></param>
        /// <param name="a"></param>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <param name="n1"></param>
        /// <param name="n2"></param>
        static void FormANkXxY(double e1, double e2, double* a, double* v1,
                               double* v2, int n1, int n2)
        {
            Parallel.For(0, n1, i =>
            {
                a[i] = e1 * v1[i];
            });

            Parallel.For(n1, n1 + n2, i =>
            {
                a[i] = e2 * v2[i - n1];
            });

        }
        #endregion

        #region FormDkXxY
        /// <summary>
        /// 
        /// </summary>
        /// <param name="D"></param>
        /// <param name="l"></param>
        /// <param name="n"></param>
        /// <param name="d"></param>
        static void FormDkXxY(double* D, double* l, int n, double d)
        {
            Parallel.For(0, n - 1, i =>
            {
                D[i] = l[i];
            });

            D[n - 1] = d;
        }
        #endregion

        #region EigenValueNew
        /// <summary>
        /// 
        /// </summary>
        /// <param name="D"></param>
        /// <param name="A"></param>
        /// <param name="n"></param>
        /// <param name="lambda"></param>
        /// <param name="epsilon"></param>
        static void EigenValueNew(double* D, double* A, int n, double* lambda, double epsilon)
        {
            //Находим собственные значения на внутренних отрезках.                        
            Parallel.For(0, n - 2, i =>
            {
                lambda[i + 1] = EigenValueInterval(D, A, i, n, epsilon);
            });

            //Находим собственные значения на граничных интервалах.     
            EigenValueInfinit(D, A, lambda, n - 1);
        }
        #endregion

        #region EigenVectorsNew
        /// <summary>
        /// 
        /// </summary>
        /// <param name="V"></param>
        /// <param name="lambda"></param>
        /// <param name="E"></param>
        /// <param name="D"></param>
        /// <param name="n"></param>
        /// <param name="up"></param>
        static void EigenVectorsNew(double* V, double* lambda, double* E, double* D, int n, bool up)
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
                vector[1, j] = ((BigDouble)(lambda[j] - *D)) * vector[0, j] / (BigDouble)(*E);


                double* d = D + 1, e = E;
                fixed (BigDouble* v = &vector[2, j])
                    for (BigDouble* st = v, sp = v + (n - 3) * n + 1; st < sp; st += n, d++, e++)
                    {
                        *st = ((BigDouble)(lambda[j] - *d) * *(st - n) - (BigDouble)(*e) * *(st - n - n)) / (BigDouble)(*(e + 1));
                    }
            });

            //Нормируем вектора
            fixed (BigDouble* v = &vector[0, 0])
                NormVectorMatrixBDnew(v, n, V, up);
        }
        #endregion

        #region NormVectorMatrixBDnew
        /// <summary>
        /// 
        /// </summary>
        /// <param name="d"></param>
        /// <param name="n"></param>
        /// <param name="vector"></param>
        /// <param name="up"></param>
        static void NormVectorMatrixBDnew(BigDouble* d, int n, double* vector, bool up)
        {
            Parallel.For(0, n, j =>
            {
                BigDouble temp = new BigDouble(0.0, 1);

                for (BigDouble* st = d + j, sp = st + (n - 1) * n + 1; st < sp; st += n)
                    temp += BigDouble.Pow(*st, 2);

                temp = BigDouble.Sqrt(temp);

                if (up)
                    vector[j] = (double)(*(d + j) / temp);
                else
                    vector[j] = (double)(*(d + j + (n - 1) * n) / temp);
            });
        }
        #endregion

        #region FormANkAddNew
        /// <summary>
        /// 
        /// </summary>
        /// <param name="AN"></param>
        /// <param name="D"></param>
        /// <param name="V"></param>
        /// <param name="L"></param>
        /// <param name="N"></param>
        /// <param name="e"></param>
        /// <param name="d"></param>
        static void FormANkAddNew(double* AN, double* D, double* V, double* L, int N, double e, double d)
        {
            Parallel.For(0, N - 1, i =>
            {
                AN[i] = e * V[i];
                D[i] = L[i];
            });

            D[N - 1] = d;
        }
        #endregion

        #region EigenValueAndVectorAddNew
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="up"></param>
        /// <param name="lambda"></param>
        /// <param name="l_old"></param>
        /// <param name="vector"></param>
        /// <param name="v_old"></param>
        /// <param name="e"></param>
        /// <param name="d"></param>
        /// <param name="epsilon"></param>
        static void EigenValueAndVectorAddNew(int n, bool up, double* lambda, double* l_old, double* vector, 
                                               double* v_old, double e, double d, double epsilon)
        {
            //Выделяем последний (новый) столбец.
            var AN = new double[n - 1];

            var D = new double[n];

            fixed (double* _a = AN, _d = D)
                FormANkAddNew(_a, _d, v_old, l_old, n, e, d);

            fixed (double* _d = D, _a = AN)
                EigenValue(_d, _a, n, lambda, vector, epsilon, up);
        }
        #endregion


        #region smatrixevvTDnew
        /// <summary>
        /// 
        /// </summary>
        /// <param name="dA"></param>
        /// <param name="eA"></param>
        /// <param name="n"></param>
        /// <param name="lambda"></param>
        /// <param name="vector"></param>
        /// <param name="epsilon"></param>
        /// <param name="vectorsFound"></param>
        /// <returns></returns>
        public static bool smatrixevvTDnew(double[] d, double[] e, int n, ref double[] lambda,
                                        ref double[,] vector, double epsilon, bool vectorsFound)
        {
            //Текущий размер вычисляемых матриц
            int del = 2;
            //Количество матриц текущего размера
            int ni = (n + 1) / (del + 1);
            //Количество оставшихся элементов
            int tale = n + 1 - ni * (del + 1);
            //Массив верхних или нижних строк собственых
            //векторов матриц текущего размера
            double[][] vectors = new double[ni][];
            double[][] lambdas = new double[ni][];

            //Вычисляем собственные числа и вектора
            //всех матриц, размерности 2
            for (int i = 0, j = 1; j <= ni; i += (del + 1), j++)
            {
                double[,] Ak = new double[del, del];
                Ak[0, 0] = d[i];
                Ak[0, 1] = e[i];
                Ak[1, 0] = e[i];
                Ak[1, 1] = d[i + 1];

                double[] L = new double[del];
                double[] V = new double[del];
                fixed (double* a = Ak, l = L, v = V)
                    EigenValueAndVector2x2TD(a, l, v, j % 2 == 0);

                lambdas[j - 1] = L;
                vectors[j - 1] = V;
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
                        EigenValueAndVectorAddNew(del + i, i == tale, l, l_old, v, v_old, e[n - tale - 2 + i],
                                                   d[n - tale - 1 + i], epsilon);

                    lambdas[ni - 1] = L;
                    if (i != tale)
                        vectors[ni - 1] = V;
                }

                fixed (double* newL = lambdas[ni - 1], newV = vectors[ni - 1], _e = &e[n - (del + tale)], _d = &d[n - (del + tale)])
                    EigenVectorsNew(newV, newL, _e, _d, del + tale, true);
            }

            //Обновляем значения
            del = del + del + 1;
            ni = (n + 1) / (del + 1);
            tale = n + 1 - ni * (del + 1);

            while (del <= n)
            {
                double[][] newVectors = new double[ni][];
                double[][] newLambdas = new double[ni][];

                int k = 0;

                for (int i = 0, j = 0; j < ni; i += (del + 1), j++, k += 2)
                {
                    //Вычисление чисел и векторов
                    //Если j % 2 = 0, записываем первую строку
                    //матрицы сосбтвенных векторов, иначе - последнюю
                    int poz = del / 2 + i;
                    newLambdas[j] = new double[del];
                    newVectors[j] = new double[del];
                    fixed (double* oldL1 = lambdas[k], oldV1 = vectors[k],
                           oldL2 = lambdas[k + 1], oldV2 = vectors[k + 1],
                           newL = newLambdas[j])
                        EigenValueAndVectorXxX(d[poz], e[poz - 1], e[poz], i, del, oldL1, oldV1, oldL2,
                                               oldV2, newL, epsilon);
                    fixed (double* newL = newLambdas[j], newV = newVectors[j], _e = &e[i], _d = &d[i])
                        EigenVectorsNew(newV, newL, _e, _d, del, (j + 1) % 2 == 0);

                }

                //Если количество матриц не четное,
                //то нужно совместить последнюю матрицу
                //размерности del с остатком
                if (ni % 2 != 0 && tale > 0)
                {
                    if (tale < 3)
                    {
                        //Посчитать остатток методом Add
                        for (int i = 1; i <= tale; i++)
                        {
                            double[] L = new double[del + i];
                            double[] V = new double[del + i];
                            fixed (double* l = L, v = V, l_old = newLambdas[ni - 1], v_old = newVectors[ni - 1])
                                EigenValueAndVectorAddNew(del + i, i == tale, l, l_old, v, v_old, e[n - tale - 2 + i],
                                                           d[n - tale - 1 + i], epsilon);

                            newLambdas[ni - 1] = L;
                            if (i != tale)
                                newVectors[ni - 1] = V;
                        }
                    }
                    else
                    {
                        double[] tempL = new double[del + tale];
                        double[] tempV = new double[del + tale];
                        //Посчитать как 
                        fixed (double* oldL = lambdas[lambdas.Length - 1], oldV = vectors[vectors.Length - 1],
                               newL = newLambdas[ni - 1], newV = newVectors[ni - 1],
                               L = tempL)
                            EigenValueAndVectorXxY(d[n - tale], e[n - tale - 1], e[n - tale], del, tale - 1,
                                                   newL, newV, oldL, oldV, L, epsilon);

                        newLambdas[ni - 1] = tempL;

                    }

                    fixed (double* newL = newLambdas[ni - 1], newV = newVectors[ni - 1], _e = &e[n - (del + tale)], _d = &d[n - (del + tale)])
                        EigenVectorsNew(newV, newL, _e, _d, del + tale, true);
                }

                lambdas = newLambdas;
                vectors = newVectors;

                //Обновляем значения
                del = del + del + 1;
                ni = (n + 1) / (del + 1);
                tale = n + 1 - ni * (del + 1);
            }

            lambda = lambdas[0];

            if (vectorsFound)
            {
                //Создаем матрицу собственных векторов.
                vector = new double[n, n];
                //Вычисляем собственные вектора
                //трехдиагональной матрицы.
                fixed (double* _d = d, _e = e, l = lambda, v = vector)
                    EigenVectors(_d, _e, l, n, v);
            }

            return true;
        }
        #endregion
    }
}

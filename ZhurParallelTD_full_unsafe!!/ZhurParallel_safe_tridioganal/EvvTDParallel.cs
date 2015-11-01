using System;
using System.Threading;
using System.Threading.Tasks;

namespace ZhurParallelTDusF
{
    public unsafe partial class EigenValueAndVectorProblem
    {
        #region EigenValueAndVectorXxY
        /// <summary>
        /// Метод вычисляет собственные числа ограниченно-диагональной матрицы
        /// образуемой из двух матриц ("верхней" и "нижней") с предыдущего
        /// шага вычисления различной размерности.
        /// </summary>
        /// <param name="d">Срдений элемент ограниченно-диагональной матрицы.</param>
        /// <param name="e1">Верхний побочный элемент ограниченно-диагональной 
        /// матрицы.</param>
        /// <param name="e2">Нижний побочный элемент ограниченно-диагональной
        /// матрицы.</param>
        /// <param name="n1">Размерность верхней матрицы.</param>
        /// <param name="n2">Размерность нижней матрицы.</param>
        /// <param name="oldL1">Указатель на массив, составляющий собственные
        /// числа верхней матрицы.</param>
        /// <param name="oldV1">Указатель на массив, составляющий нижнюю строку
        /// матрицы собственных вектора верхней матрицы.</param>
        /// <param name="oldL2">Указатель на массив, составляющий собственные
        /// числа нижней матрицы.</param>
        /// <param name="oldV2">Указатель на массив, составляющий верхнюю строку
        /// матрицы собственных вектора верхней матрицы.</param>
        /// <param name="newL">Указатель на массив, составляющий собственные
        /// числа ограниченно-диагональной матрицы.</param>
        /// <param name="epsilon">Точность вычисления.</param>
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
        /// Метод вычисляет собственные числа ограниченно-диагональной матрицы
        /// образуемой из двух матриц ("верхней" и "нижней") с предыдущего
        /// шага вычисления одной размерности.
        /// </summary>
        /// <param name="d">Срдений элемент ограниченно-диагональной матрицы.</param>
        /// <param name="e1">Верхний побочный элемент ограниченно-диагональной 
        /// матрицы.</param>
        /// <param name="e2">Нижний побочный элемент ограниченно-диагональной
        /// матрицы.</param>
        /// <param name="n">Размерность ограниченно-диагональной матрицы.</param>
        /// <param name="oldL1">Указатель на массив, составляющий собственные
        /// числа верхней матрицы.</param>
        /// <param name="oldV1">Указатель на массив, составляющий нижнюю строку
        /// матрицы собственных вектора верхней матрицы.</param>
        /// <param name="oldL2">Указатель на массив, составляющий собственные
        /// числа нижней матрицы.</param>
        /// <param name="oldV2">Указатель на массив, составляющий верхнюю строку
        /// матрицы собственных вектора верхней матрицы.</param>
        /// <param name="newL">Указатель на массив, составляющий собственные
        /// числа ограниченно-диагональной матрицы.</param>
        /// <param name="epsilon">Точность вычисления.</param>
        static void EigenValueAndVectorXxX(double d, double e1, double e2, int n, double* oldL1, double* oldV1,
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
        /// Метод формирует массив собственных чисел из массивов собственных 
        /// чисел верхней и нижней матриц.
        /// </summary>
        /// <param name="L">Указатель на массив, содержащий составленные
        /// собственные числа.</param>
        /// <param name="L1">Указатель на массив, содержащий собственные числа
        /// верхней матрицы.</param>
        /// <param name="L2">Указатель на массив, содержащий собственные числа
        /// нижней матрицы.</param>
        /// <param name="n1">Размерность верхней матрицы.</param>
        /// <param name="n2">Размерность нижней матриццы.</param>
        static void FormL(double* L, double* L1, double* L2, int n1, int n2)
        {
            Parallel.For(0, n1, i =>
            //for(int i = 0; i < n1; i++)
            {
                L[i] = L1[i];
            });

            Parallel.For(n1, n1 + n2, i =>
            //for (int i = n1; i < n1 + n2; i++)
            {
                L[i] = L2[i - n1];
            });
        }
        #endregion

        #region FormANkXxY
        /// <summary>
        /// Метод формирует столбец ограниченно-диагональной матрицы.
        /// </summary>
        /// <param name="e1">Верхний внедиагональный элемент 
        /// ограниченно-диагональной матрицы.</param>
        /// <param name="e2">Нижний внедиагональный элемент
        /// ограниченно-диагональной матрицы.</param>
        /// <param name="a">Указатель на массив, составляющий
        /// столбец ограниченно-диагональной матрицы.</param>
        /// <param name="v1">Указатель на массив, составляющий последнюю
        /// строку матрицы собственных векторов, находящейся сверху.</param>
        /// <param name="v2">Указатель на массив, составляющий первую
        /// строку матрицы собственных векторов, находящейся внизу.</param>
        /// <param name="n1">Размерность верхней матрицы.</param>
        /// <param name="n2">Размерность нижней матрицы.</param>
        static void FormANkXxY(double e1, double e2, double* a, double* v1,
                               double* v2, int n1, int n2)
        {
            Parallel.For(0, n1, i =>
            //for(int i = 0; i < n1; i++)
            {
                a[i] = e1 * v1[i];
            });

            Parallel.For(n1, n1 + n2, i =>
            //for (int i = n1; i < n1 + n2; i++)
            {
                a[i] = e2 * v2[i - n1];
            });

        }
        #endregion

        #region FormDkXxY
        /// <summary>
        /// Метод формирует главную диагональ ограниченно-диагональной матрицы.
        /// </summary>
        /// <param name="D">Указатель на массив, составляющий главную диагональ
        /// ограниченно-диагональной матрицы.</param>
        /// <param name="l">Указатель на массив, составляющий собственные числа
        /// из предыдущего шага.</param>
        /// <param name="n">Размерность ограниченно-диагональной матрицы.</param>
        /// <param name="d">Последний диагональный элемент ограниченно-диагональной 
        /// матрицы.</param>
        static void FormDkXxY(double* D, double* l, int n, double d)
        {
            Parallel.For(0, n - 1, i =>
            //for (int i = 0; i < n - 1; i++)
            {
                D[i] = l[i];
            });

            D[n - 1] = d;
        }
        #endregion

        #region EigenValueNew
        /// <summary>
        /// Метод вычисляет собственные числа ограниченно-диагональной матрицы.
        /// </summary>
        /// <param name="D">Указатель на массив, составляющий главную диагональ 
        /// ограниченно-диагональной матрицы.</param>
        /// <param name="A">Указатель на массив, составляющий столбец 
        /// ограниченно-диагональной матрицы.</param>
        /// <param name="n">Размерность ограниченно-диагональной матрицы.</param>
        /// <param name="lambda">Указатель на массив, соствляющий полученные
        /// собственные числа.</param>
        /// <param name="epsilon">Точность вычисления.</param>
        static void EigenValueNew(double* D, double* A, int n, double* lambda, double epsilon)
        {
            //Находим собственные значения на внутренних отрезках.                        
            Parallel.For(0, n - 2, i =>
            //for (int i = 0; i < n - 2; i++)
            {
                lambda[i + 1] = EigenValueInterval(D, A, i, n, epsilon);
            });

            //Находим собственные значения на граничных интервалах.     
            EigenValueInfinit(D, A, lambda, n - 1);
        }
        #endregion

        #region EigenVectorsNew
        /// <summary>
        /// Метод вычисляет собственные вектора трехдиагональной матрицы и возвращает выбранную строку.
        /// </summary>
        /// <param name="V">Указатель на массив, содержащий строку матрицы собственных в зависимости
        ///  от <paramref name="up"/>.</param>
        /// <param name="lambda">Указатель на массив, составляющий собственные числа трехдиагональной матрицы.</param>
        /// <param name="E">Указатель на массив, составляющий побочную диагональ трехдиагональной матрицы.</param>
        /// <param name="D">Указатель на массив, составляющий главную диагональ трехдиагональной матрицы.</param>
        /// <param name="n">Размерность трехдиагональной матрицы.</param>
        /// <param name="up">Параметр, указывающий какую из строк (верхнюю или нижнюю) необходимо вернуть.</param>
        static void EigenVectorsNew(double* V, double* lambda, double* E, double* D, int n, bool up)
        {
            //По скольку, при вычислении собственных векторов
            //трехдиагональной матрицы, значения временно могут
            //выходить за рамки значений типа Double,
            //создаем матрицу специального типа.            
            var tVector = new double[n * n];

            //Вычисляем собственные вектора
            //трехдиагональной матрицы с
            //помощью метода прогонки
            var epsi = new double[n * n];
            var etha = new double[n * n];
            fixed (double* tV = tVector, _epsi = epsi, _etha = etha)
            {
                var eps = _epsi;
                var et = _etha;

                var vector = tV;

                Parallel.For(0, n, j =>
                //for (int j = 0; j < n; j++)
                {
                    var d = Math.Sqrt(1.0 / n);

                    var ep = eps + n * j;
                    var eta = et + n * j;

                    ep[1] = -(*E / (*D - lambda[j]));
                    eta[1] = d / (*D - lambda[j]);

                    //for (int i = 1; i < n - 1; i++) 
                    for (double* step = ep + 2, steta = eta + 2, en = ep + n, e = E, _d = D + 1, l = lambda + j; step < en; steta++, step++, e++, _d++)
                    {
                        //ep[i + 1] = -E[i] / (D[i] - lambda[j] + E[i - 1] * ep[i]);
                        *step = -(*(e + 1) / (*_d - *l + (*e * *(step - 1))));
                        //eta[i + 1] = (d - E[i - 1] * eta[i]) / (D[i] - lambda[j] + E[i - 1] * ep[i]);
                        *steta = (d - (*e * *(steta - 1))) / (*_d - *l + (*e * *(step - 1)));
                    }

                    var lastEta = (d - E[n - 2] * eta[n - 1]) / (D[n - 1] - (lambda[j] - 1e-6) + E[n - 2] * ep[n - 1]);
                    vector[(n - 1) * n + j] = lastEta;
                    var norm = Sqr(lastEta);

                    //for (int i = n - 1; i > 1; i--)
                    for (double* v = vector + (n - 1) * n + j, _ep = ep + n - 1, _eta = eta + n - 1, en = vector + j; v > en; v -= n, _ep--, _eta--)
                    {
                        *(v - n) = *_ep * *v + *_eta;
                        //norm += Sqr(*(v - n));
                    }

                    //norm = Math.Sqrt(norm);

                    ////for (int i = 0; i < n; i++)
                    //for (double* v = vector + j, en = vector + n * n; v < en; v += n)
                    //    *v = *v / norm;
                });

                //    //Второй проход для большей точности
                //    Parallel.For(0, n, j =>
                //    //for (int j = 0; j < n; j++)
                //    {
                //        var ep = eps + n * j;
                //        var eta = et + n * j;

                //        ep[1] = -(*E / (*D - lambda[j]));
                //        eta[1] = vector[j] / (*D - lambda[j]);

                //        //for (int i = 1; i < n - 1; i++)
                //        for (double* step = ep + 2, steta = eta + 2, en = ep + n, e = E, _d = D + 1, l = lambda + j, v = vector + n + j; step < en; steta++, step++, e++, _d++, v += n)
                //        {
                //            //ep[i + 1] = E[i] / (D[i] - lambda[j] + E[i - 1] * ep[i]);
                //            *step = -(*(e + 1) / (*_d - *l + (*e * *(step - 1))));
                //            //eta[i + 1] = (vector[i * n + j] - E[i - 1] * eta[i]) / (D[i] - lambda[j] + E[i - 1] * ep[i]);
                //            *steta = (*v - (*e * *(steta - 1))) / (*_d - *l + (*e * *(step - 1)));
                //        }

                //        vector[(n - 1) * n + j] = (vector[(n - 1) * n + j] - E[n - 2] * eta[n - 1]) / (D[n - 1] - (lambda[j] - 1e-6) + E[n - 2] * ep[n - 1]);

                //        for (double* v = vector + (n - 1) * n + j, _ep = ep + n - 1, _eta = eta + n - 1, en = vector + j; v > en; v -= n, _ep--, _eta--)
                //        {
                //            *(v - n) = *_ep * *v + *_eta;
                //        }

                //    });
            }

            epsi = null;
            etha = null;

            //Нормируем вектора
            //fixed (BigDouble* v = &vector[0, 0])
            //    NormVectorMatrixBDnew(v, n, V, up);
            fixed (double* vector = tVector)
                NormVector(vector, n, V, up);

            tVector = null;
        }

        /// <summary>
        /// Метод нормирует матрицу собственных векторов и возвращает выбранную строку результата.
        /// </summary>
        /// <param name="vector">Указатель на матрицу собственных векторов.</param>
        /// <param name="n">Размерность матрицы.</param>
        /// <param name="V">Указатель на массив, составляющий в зависимости от <paramref name="up"/>
        /// троку нормированной матрицы векторов. </param>
        /// <param name="up">Параметр, указывающий на то, какую строку (верхнюю или нижнюю) нужно вернуть.</param>
        static void NormVector(double* vector, int n, double* V, bool up)
        {
            int add;
            if (up)
                add = 0;
            else
                add = n * (n - 1);

            Parallel.For(0, n, j =>
                {
                    var temp = 0.0d;

                    for (double* st = vector + j, sp = vector + n * n; st < sp; st += n)
                        temp += Sqr(*st);

                    temp = Math.Sqrt(temp);
                                        
                    V[j] = *(vector + add + j) / temp;
                    
                });
        }
        #endregion
        
        #region FormANkAddNew
        /// <summary>
        /// Метод формирует столбец ограниченно-диагональной матрицы.
        /// </summary>
        /// <param name="AN">Указатель на массив, составляющий столбец ограниченно-диагональной матрицы.</param>
        /// <param name="D">Указатель на массив, составляющий главную диагональ ограниченно-диагональной матрицы.</param>
        /// <param name="V">Указатель на массив, составляющий строку собственных векторов из предыдущего шага.</param>
        /// <param name="L">Указатель на массив, составляющий собственные числа из предыдущего шага.</param>
        /// <param name="N">Размерность ограниченно-диагональной матрицы.</param>
        /// <param name="e">Внедиагональный элемент ограниченно-диагональной матрицы.</param>
        /// <param name="d">Последний элемент ограниченно-диагональной матрицы.</param>
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
        /// Метод вычисляет собственные числа и вектора дополненной матрицы.
        /// </summary>
        /// <param name="n">Размерность матрицы.</param>
        /// <param name="up">Параметр, указывающий на то какую строку собственных векторов
        /// необходимо вернуть.</param>
        /// <param name="lambda">Указатель на массив полученных собственных чисел.</param>
        /// <param name="l_old"> Указатель на массив собственных чисел на предыдущем шаге.</param>
        /// <param name="vector">В зависимости от <paramref name="up"/> возвращает указатель 
        /// на массив, составляющий верхнюю или нижнюю строку полученной матрицы собственных векторов.</param>
        /// <param name="v_old"></param>
        /// <param name="e">Указатель на массив составляющий побочную диагональ дополненной матрицы.</param>
        /// <param name="d">Указатель на массив составляющий главную диагональ дополненной матрицы.</param>
        /// <param name="epsilon">Точность вычисления.</param>
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
        /// Основной метод вычисления собственных чисел и векторов трехдиагональной матрицы.
        /// </summary>
        /// <param name="d">Массив, составляющий главную диагональ трехдиагональной матрицы.</param>
        /// <param name="e">Массив, составляющий побочную диагональ трехдиагональной матрицы.</param>
        /// <param name="n">Размерность матрицы.</param>
        /// <param name="lambda">Массив полученных собственных чисел.</param>
        /// <param name="vector">Матрица полученных собственных векторов (Если <paramref name="vectorsFound"/> = false).</param>
        /// <param name="epsilon">Точность вычисления.</param>
        /// <param name="vectorsFound">Параметр, указывающий на необходимость вычисления собственных векторов.</param>
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
            //for (int j = 0; j < ni; j++)
            Parallel.For(0, ni, j =>
            {
                var i = j * (del + 1);
                double[] Ak = new double[del * del];
                Ak[0] = d[i];
                Ak[1] = e[i];
                Ak[2] = e[i];
                Ak[3] = d[i + 1];

                double[] L = new double[del];
                double[] V = new double[del];
                fixed (double* a = Ak, l = L, v = V)
                    EigenValueAndVector2x2TD(a, l, v, (j + 1) % 2 == 0);

                lambdas[j] = L;
                vectors[j] = V;

                L = null;
                V = null;
            });

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
                    L = null;
                    if (i != tale)
                        vectors[ni - 1] = V;
                }

                vectors[ni - 1] = new double[del + tale];
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

                //for (int j = 0; j < ni; j++)
                Parallel.For(0, ni, j =>
                {

                    var i = j * (del + 1);
                    var k = j * 2;
                    //Вычисление чисел и векторов
                    //Если j % 2 = 0, записываем первую строку
                    //матрицы сосбтвенных векторов, иначе - последнюю
                    var poz = del / 2 + i;
                    newLambdas[j] = new double[del];
                    newVectors[j] = new double[del];
                    fixed (double* oldL1 = lambdas[k], oldV1 = vectors[k],
                           oldL2 = lambdas[k + 1], oldV2 = vectors[k + 1],
                           newL = newLambdas[j])
                        EigenValueAndVectorXxX(d[poz], e[poz - 1], e[poz], del, oldL1, oldV1, oldL2,
                                               oldV2, newL, epsilon);
                    if (del != n)
                        fixed (double* newL = newLambdas[j], newV = newVectors[j], _e = &e[i], _d = &d[i])
                            EigenVectorsNew(newV, newL, _e, _d, del, (j + 1) % 2 == 0);

                });

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
                            L = null;
                            if (i != tale)
                                newVectors[ni - 1] = V;
                        }
                    }
                    else
                    {
                        double[] tempL = new double[del + tale];                        
                        fixed (double* oldL = lambdas[lambdas.Length - 1], oldV = vectors[vectors.Length - 1],
                               newL = newLambdas[ni - 1], newV = newVectors[ni - 1],
                               L = tempL)
                            EigenValueAndVectorXxY(d[n - tale], e[n - tale - 1], e[n - tale], del, tale - 1,
                                                   newL, newV, oldL, oldV, L, epsilon);

                        newLambdas[ni - 1] = tempL;
                        tempL = null;

                    }

                    if (del + tale != n)
                    {
                        newVectors[ni - 1] = new double[del + tale];
                        fixed (double* newL = newLambdas[ni - 1], newV = newVectors[ni - 1], _e = &e[n - (del + tale)], _d = &d[n - (del + tale)])
                            EigenVectorsNew(newV, newL, _e, _d, del + tale, true);
                    }
                }
                else
                {
                    if (tale > 0)
                    {
                        double[][] tempVectors = new double[ni + 1][];
                        double[][] tempLambdas = new double[ni + 1][];

                        Parallel.For(0, ni, i =>
                        {
                            tempLambdas[i] = newLambdas[i];
                            tempVectors[i] = newVectors[i];
                        });

                        tempVectors[ni] = vectors[vectors.Length - 1];
                        tempLambdas[ni] = lambdas[lambdas.Length - 1];

                        newLambdas = tempLambdas;
                        newVectors = tempVectors;
                    }

                }

                lambdas = newLambdas;
                vectors = newVectors;

                newLambdas = null;
                newVectors = null;
                
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

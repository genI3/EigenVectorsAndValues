﻿using System;
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
                    if (tt == 0.0d)
                        t = 0.0d;
                    else
                        t = Sqr(*_a / tt);

                    temp += t;
                }

                if (up)
                {
                    var del = lambda[j] - *D;
                    if (del == 0.0d)
                        V[j] = 0.0d;
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
            //Вычисляем собственные вектора
            //трехдиагональной матрицы в BigDouble

            var epsi = new double[n * n];
            var etha = new double[n * n];

            fixed (double* _epsi = epsi, _etha = etha)
            {
                var eps = _epsi;
                var et = _etha;

                Parallel.For(0, n, j =>
                {
                    var d = Math.Sqrt(1.0 / n);

                    var ep = eps + n * j;
                    var eta = et + n * j;

                    ep[1] = -(*E / (*D - L[j]));
                    eta[1] = d / (*D - L[j]);

                    //for (int i = 1; i < n - 1; i++)
                    for (double* step = ep + 2, steta = eta + 2, en = ep + n, e = E, _d = D + 1, l = L + j; step < en; steta++, step++, e++, _d++)
                    {
                        //ep[i + 1] = -E[i] / (D[i] - L[j] + E[i - 1] * ep[i]);
                        *step = -(*(e + 1) / (*_d - *l + (*e * *(step - 1))));
                        //eta[i + 1] = (d - E[i - 1] * eta[i]) / (D[i] - L[j] + E[i - 1] * ep[i]);
                        *steta = (d - (*e * *(steta - 1))) / (*_d - *l + (*e * *(step - 1)));
                    }

                    //Получаем вектор
                    res[(n - 1) * n + j] = (d - E[n - 2] * eta[n - 1]) / (D[n - 1] - (L[j] - 1e-6) + E[n - 2] * ep[n - 1]);
                    var norm = Sqr(res[(n - 1) * n + j]);

                    for (double* v = res + (n - 1) * n + j, _ep = ep + n - 1, _eta = eta + n - 1, en = res + j; v > en; v -= n, _ep--, _eta--)
                    {
                        *(v - n) = *_ep * *v + *_eta;
                        norm += Sqr(*(v - n));
                    }

                    norm = Math.Sqrt(norm);

                    //Нормируем его
                    for (double* v = res + j, en = res + n * n; v < en; v += n)
                        *v = *v / norm;

                });

                //Второй проход для большей точности
                Parallel.For(0, n, j =>
                //for (int j = 0; j < n; j++)
                {
                    var ep = eps + n * j;
                    var eta = et + n * j;

                    ep[1] = -(*E / (*D - L[j]));
                    eta[1] = res[j] / (*D - L[j]);

                    //for (int i = 1; i < n - 1; i++)
                    for (double* step = ep + 2, steta = eta + 2, en = ep + n, e = E, _d = D + 1, l = L + j, v = res + n + j; step < en; steta++, step++, e++, _d++, v += n)
                    {
                        //ep[i + 1] = -E[i] / (D[i] - L[j] + E[i - 1] * ep[i]);
                        *step = -(*(e + 1) / (*_d - *l + (*e * *(step - 1))));
                        //eta[i + 1] = (res[i * n + j] - E[i - 1] * eta[i]) / (D[i] - L[j] + E[i - 1] * ep[i]);
                        *steta = (*v - (*e * *(steta - 1))) / (*_d - *l + (*e * *(step - 1)));
                    }

                    res[(n - 1) * n + j] = (E[n - 2] * eta[n - 1] - res[(n - 1) * n + j]) / (D[n - 1] - (L[j] - 1e-6) + E[n - 2] * ep[n - 1]);
                    var norm = Sqr(res[(n - 1) * n + j]);

                    for (double* v = res + (n - 1) * n + j, _ep = ep + n - 1, _eta = eta + n - 1, en = res + j; v > en; v -= n, _ep--, _eta--)
                    {
                        *(v - n) = *_ep * *v + *_eta;
                        norm += Sqr(*(v - n));
                    }

                    norm = Math.Sqrt(norm);

                    for (double* v = res + j, en = res + n * n; v < en; v += n)
                        *v = *v / norm;
                });
            }

            etha = null;
            epsi = null;
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

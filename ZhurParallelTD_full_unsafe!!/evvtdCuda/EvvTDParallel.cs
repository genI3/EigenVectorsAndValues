using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Cudafy;
using Cudafy.Host;
using Cudafy.Translator;

namespace ZhurParallelTDCuda
{
    public class evvtdcuda
    {
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
            //CUDA initialize
            CudafyModule km = CudafyTranslator.Cudafy();
            GPGPU gpu = CudafyHost.GetDevice(CudafyModes.Target, CudafyModes.DeviceId);
            gpu.LoadModule(km);

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

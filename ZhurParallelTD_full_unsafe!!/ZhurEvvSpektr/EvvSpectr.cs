using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZhurEvvSpektr
{
    public static class evvspectr
    {
        public static void evvspectrc(double[,] X, double[] x, double[,] R, double[] l_R,
                                      double[,] v_R, out double[] l, out double[,] v, int N, int M)
        {
            var mid_XN = new double[M];
            var mid_XN1 = new double[M];

            for (int i = 0; i < M; i++)
            {
                var tmp = 0.0d;

                for (int j = 0; j < N; j++)
                    tmp += X[i, j];

                mid_XN[j] = tmp / N;
                mid_XN1[j] = (tmp + x[i]) / (N + 1);
            }

            var sigmaN  = new double[M];
            var sigmaN1 = new double[M];
            
            for (int i = 0; i < M; i++)
            {
                var tmp = 0.0d;

                for (int j = 0; j < N; j++)
                    tmp += Sqr(X[i, j] - mid_XN[i]);

                sigmaN[i]  = Math.Sqrt(tmp / (N - 1));
                sigmaN1[i] = Math.Sqrt((tmp + Sqr(x[i] - mid_XN1[i])) / N);
            }

            var ra = new double[M];
            var ya = new double[M];

            for (int i = 0; i < M; i++)
            {
                ra[i] = (Math.Sqrt(1.0d + Sqr(N)) / (1 + N)) *
                        (x[i] - mid_XN[i]);
                ya[i] = ra[i] * sigmaN[i];
            }

            var La = new double[M];
            var cN = (N - 1) * 1.0d / N;

            for (int i = 0; i < M; i++)
            {
                La[i] = cN * l_R[i];
            }

            EigenValuesAndVectors(ya, La, out l, out v, M);

            for (int i = 0; i < M; i++)
                l[i] = l[i] * Sqr(sigmaN[i] / sigmaN1[i]);

            EgenVectors(sigmaN1, sigmaN, N, v_R, v);

        }

        /// <summary>
        /// [S(N+1)]^-1 * S(N) * v_R * v
        /// </summary>
        /// <param name="sigmaN1"></param>
        /// <param name="sigmaN"></param>
        /// <param name="n"></param>
        /// <param name="v_R"></param>
        /// <param name="v"></param>
        static void EgenVectors(double[] sigmaN1, double[] sigmaN, int n, double[,] v_R, double[,] v)
        {
            throw new NotImplementedException();
        }


        /// <summary>
        /// Solve equation
        /// </summary>
        /// <param name="ya"></param>
        /// <param name="la"></param>
        /// <param name="l"></param>
        /// <param name="v"></param>
        /// <param name="m"></param>
        static void EigenValuesAndVectors(double[] ya, double[] la, out double[] l, out double[,] v, int m)
        {
            throw new NotImplementedException();
        }

        static double Sqr(double m)
        {
            return m * m;
        }
     
        static int Sqr(int m)
        {
            return m * m;
        }       
    }
}

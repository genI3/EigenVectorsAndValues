using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ZhurParallelTDusF;

namespace NewTest
{
    class Program
    {
        static void Main(string[] args)
        {
            int n = 15000;
            
            double[] d;
            double[] e;

            double[] lambda = new double[2];
            double[,] vector = new double[2, 2];

            //RandomTridiagonalMatrix(n, -13.6d, 10.3d, 0);


            //using (StreamReader sr = new StreamReader("arrays.txt"))
            //{
            //    //BinaryReader br = new BinaryReader(fs);
            //    var data = sr.ReadToEnd().Split(new char[] { '\n' });
            //    int n1 = int.Parse(data[0]);// ReadInt32();

            //    d = new double[n1];
            //    for (int i = 0; i < n1; i++)
            //        d[i] = double.Parse(data[i + 1]);//br.ReadDouble();

            //    e = new double[n1 - 1];
            //    for (int i = 0; i < n1 - 1; i++)
            //        e[i] = double.Parse(data[i + 1 + n1]);
            //}

            RandomTridiagonalMatrix(n, -13.6d, 10.3d, 0, out d, out e);

            Console.WriteLine("Исходная трехдиагональная матрица размерности {0}", n);
            Console.Write("Главная диагональ: ");
            for (int i = 0; i < 6; i++)
                Console.Write("{0:0.000} ", d[i]);
            Console.WriteLine();

            Console.Write("Побочная диагональ: ");
            for (int i = 0; i < 6; i++)
                Console.Write("{0:0.000} ", e[i]);
            Console.WriteLine();

            //Console.WriteLine("Производим вычисления старым способом.");

            ///*ZhurParallelTDusF.EPlTDusF.smatrixevvTD*/
            //EigenValueAndVectorProblem.smatrixevvTD((double[])d.Clone(), (double[])e.Clone(), n, ref lambda, ref vector, 1.0e-10d, false);
            
            //Console.Write("Полученные собственные числа: ");
            //for (int i = 0; i < 10; i++)
            //    Console.Write("{0:0.00000} ", lambda[i]);


            Console.WriteLine("\nПроизводим вычисления новым способом.");

            lambda = new double[2];
            vector = new double[2, 2];

            //evvTdTemp.evvTd.smatrixTDtemp((double[])d.Clone(), (double[])e.Clone(), n, ref lambda, ref vector, 1e-9d, false);
            EigenValueAndVectorProblem.smatrixevvTDnew((double[])d.Clone(), (double[])e.Clone(), n, ref lambda, ref vector, 1e-10d, false);

            Console.Write("Полученные собственные числа: ");
            for (int i = 0; i < 10; i++)
                Console.Write("{0:0.00000} ", lambda[i]);

            Console.WriteLine("\nПроизводим вычисления через AlgLib.");

            lambda = new double[2];
            vector = new double[2, 2];

            Eigen.Evd.smatrixtdevd(ref d, e, n, 0, ref vector);

            Console.Write("Полученные собственные числа: ");
            for (int i = 0; i < 10; i++)
                Console.Write("{0:0.00000} ", d[i]);


            Console.ReadKey();

        }

        public static void RandomTridiagonalMatrix(int n, double minVal, double maxVal, int seed)
        {
            double[] d = new double[n];
            double[] e = new double[n - 1];

            Random ran = new Random(seed);
            Parallel.For(0, n - 1, i =>
            {
                d[i] = (maxVal - minVal) * ran.NextDouble() + minVal;
                e[i] = (maxVal - minVal) * ran.NextDouble() + minVal;
             
            });

            d[n - 1] = (maxVal - minVal) * ran.NextDouble() + minVal;

            if (File.Exists("arrays.txt"))
                File.Delete("arrays.txt");

            //FileStream fs = new FileStream("arrays.txt", FileMode.OpenOrCreate);
            using (var sr = new StreamWriter("arrays.txt"))
            {
                // BinaryWriter bw = new BinaryWriter(fs);
                //StreamWriter sr = new StreamWriter(fs);

                sr.WriteLine(n.ToString());
                for (int i = 0; i < n; i++)
                {
                    sr.WriteLine(d[i].ToString());
                }

                for (int i = 0; i < n - 1; i++)
                {
                    sr.WriteLine(e[i].ToString());
                }
            }
        }

        public static void RandomTridiagonalMatrix(int n, double minVal, double maxVal, int seed, out double[] d, out double[] e)
        {
            var td = new double[n];
            var te = new double[n - 1];

            Random ran = new Random(seed);
            Parallel.For(0, n - 1, i =>
            {
                td[i] = (maxVal - minVal) * ran.NextDouble() + minVal;
                te[i] = (maxVal - minVal) * ran.NextDouble() + minVal;

            });

            td[n - 1] = (maxVal - minVal) * ran.NextDouble() + minVal;

            d = td;
            e = te;
        }
    }
}

using System;
using System.Runtime.InteropServices;

namespace ZhurParallelTDusF
{
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

            while (Math.Abs(d1.X) > 9.99d)
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

            while (Math.Abs(d1.X) > 9.99d)
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
            while (Math.Abs(d1.X) > 9.99d)
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
            while (Math.Abs(d1.X) > 9.99d)
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

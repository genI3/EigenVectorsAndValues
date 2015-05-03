using System;

namespace ZhurParallelTDusF
{
    public unsafe class smatrix2TD
    {
        const double machineepsilon = 5E-16;
        const double maxrealnumber = 1E300;
        const double minrealnumber = 1E-300;


        /*************************************************************************
       Преобразование симметричной матрицы, которая получена ее верхней или нижней
       треугольной частью к трехдиагональной матрице, используя простое ортогональное
       преобразование: Q'*A*Q=T.

       Входные параметры:
           A       -   Матрица для преобразования.
                       Массив с элементами [0..N-1, 0..N-1].
           N       -   Размерность матрицы A.
           IsUpper -   Формат хранения. Если IsUpper = True, то матрица А получена
                       ее верхней треугольной частью, нижняя треугольная часть не 
                       используется и не меняется алгоритмом, и наоборот. Если IsUpper 
                       = False.

       Выходные параметры:
           A       -   Матрицы T и Q в компактном виде (см. ниже).
           Tau     -   Массив факторов, которые формируют матрицу H(i).
                       Массив с элементами [0..N-2].
           D       -   Главная диагональ симметричной матрицы T.
                       Массив с элементами [0..N-1].
           E       -   Побочная диагональ симметричной матрицы T.
                       Массив с элементами [0..N-2].


         Если IsUpper=True, матрица Q представляет из себя продукт элементарных
         отражений

            Q = H(n-2) . . . H(2) H(0).

         Каждый H(i) имеет форму

            H(i) = I - tau * v * v'

         где tau - вещественный скаляр, и v вещественный вектор с
         v(i+1:n-1) = 0, v(i) = 1, v(0:i-1) хранится на выходе в
         A(0:i-1,i+1), и tau в TAU(i).

         Если IsUpper=False, матрица Q представляет из себя продукт элементарных
         отражений

            Q = H(0) H(2) . . . H(n-2).

         Каждый H(i) имеет форму

            H(i) = I - tau * v * v'

         где tau - вещественный скаляр, и v вещественный вектор с
         v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) зранящихся на выходе в A(i+2:n-1,i),
         и tau в TAU(i).

         Содержание A на выходе иллюстрируется следующим примером
         с n = 5:

         Если UPLO = 'U':                       Если UPLO = 'L':

           (  d   e   v1  v2  v3 )              (  d                  )
           (      d   e   v2  v3 )              (  e   d              )
           (          d   e   v3 )              (  v0  e   d          )
           (              d   e  )              (  v0  v1  e   d      )
           (                  d  )              (  v0  v1  v2  e   d  )

         где d и e обозначают диагональные и недиагональные элементы T, и vi
         обознаюает элементы вектора определения H(i).

       *************************************************************************/
        public static void smatrixtd(ref double[,] a, int n, bool isupper, ref double[] tau,
                                     ref double[] d, ref double[] e)
        {
            double alpha = 0.0d;
            double taui = 0.0d;
            double v = 0.0d;
            double[] t = new double[0];
            double[] t2 = new double[0];
            double[] t3 = new double[0];

            tau = new double[0];
            d = new double[0];
            e = new double[0];

            if (n <= 0)
            {
                return;
            }

            t = new double[n + 1];
            t2 = new double[n + 1];
            t3 = new double[n + 1];
            d = new double[n];

            if (n > 1)
            {
                tau = new double[n - 1];
            }

            if (n > 1)
            {
                e = new double[n - 1];
            }

            fixed (double* _a = a, _t = t, _e = e, _tau = tau, _d = d)
                if (isupper)
                {

                    //
                    // Уменьшаем верхний треугольник A
                    //
                    int i = n - 1;
                    for (double* stj = _a + n - 1, spj = _a + 1, spi = _a + (n - 2) * n - 1, E = _e + n - 2, D = _d + n - 1, Tau = _tau + n - 2;
                    stj >= spj;
                    stj--, spi -= (n + 1), i--, E--, D--, Tau--)
                    {

                        //
                        // Генерируем элементарное отражение H() = E - tau * v * v'
                        //
                        if (stj >= _a + 2)
                        {
                            for (double* sti = stj, ti = _t + 2; sti <= spi; sti += n, ti++)
                            {
                                *ti = *sti;
                            }
                        }
                        _t[1] = *(spi + n);
                        generatereflection(_t, i, out taui);

                        if (stj >= _a + 2)
                        {
                            //for (i_ = 0; i_ <= i - 1; i_++)
                            for (double* sti = stj, ti = _t + 2; sti <= spi; sti += n, ti++)
                            {
                                *sti = *ti;
                                //a[i_, i + 1] = t[i_ + 2];
                            }
                        }
                        *(spi + n) = _t[1];
                        *E = *(spi + n);
                        fixed (double* _t3 = t3, _t2 = t2)
                            if (taui != 0.0d)
                            {

                                //
                                // Применяем отражение H с двух сторон к A
                                //
                                *(spi + n) = 1.0d;

                                //
                                // Рассчитываем  x := tau * A * v  storing x in TAU
                                //
                                for (double* sti = stj, ti = _t + 1; sti <= spi + n; sti += n, ti++)
                                {
                                    *ti = *sti;
                                }
                                symmetricmatrixvectormultiply(_a, isupper, 0, i - 1, t, taui, _t3, n);
                                for (double* st = _tau, sp = _tau + i - 1, t3i = _t3 + 1; st <= sp; st++, t3i++)
                                {
                                    *st = *t3i;
                                }

                                //
                                // Рассчитываем  w := x - 1/2 * tau * (x'*v) * v
                                //
                                v = 0.0d;
                                for (double* sti = stj, _taui = _tau; sti <= spi + n; sti += n, _taui++)
                                {
                                    v += *_taui * *sti;
                                }
                                alpha = -(0.5d * taui * v);
                                for (double* sti = stj, _taui = _tau; sti <= spi + n; sti += n, _taui++)
                                {
                                    *_taui += alpha * *sti;
                                }

                                //
                                // Применяем преобразование как обновление rank-2:
                                //    A := A - v * w' - w * v'
                                //
                                for (double* sti = stj, _ti = _t + 1; sti <= spi + n; sti += n, _ti++)
                                {
                                    *_ti = *sti;
                                }
                                for (double* st = _tau, _t3i = _t3 + 1, sp = _tau + i; st <= sp; st++, _t3i++)
                                {
                                    *_t3i = *st;
                                }
                                symmetricrank2update(_a, isupper, 0, i - 1, t, t3, _t2, -1.0d, n);
                                *(spi + n) = *E;
                            }
                        *D = *(spi + n + n);
                        *Tau = taui;
                    }
                    *_d = *_a;
                }
                else
                {

                    //
                    // Уменьшаем нижний треугольник A
                    //
                    int i = 0;
                    for (double* stj = _a + n, spj = _a + n * n - 2, spi = _a + (n - 1) * n, E = _e, D = _d, Tau = _tau;
                    stj <= spj;
                    stj += (n + 1), spi++, i++, E++, D++, Tau++)
                    {

                        //
                        // Генерируем элементарное отражение H = E - tau * v * v'
                        //
                        for (double* sti = stj, ti = _t + 1; sti <= spi; sti += n, ti++)
                        {
                            *ti = *sti;
                            //t[i_] = a[i_ + i, i];
                        }
                        generatereflection(_t, n - i - 1, out taui);
                        for (double* sti = stj, ti = _t + 1; sti <= spi; sti += n, ti++)
                        {
                            *sti = *ti;
                        }
                        *E = *stj;
                        fixed (double* _t3 = t3, _t2 = t2)
                            if (taui != 0.0d)
                            {

                                //
                                // Применяем отражение H с двух сторон к A
                                //
                                *stj = 1.0d;

                                //
                                // Рассчитываем  x := tau * A * v  storing y in TAU
                                //
                                for (double* sti = stj, ti = _t + 1; sti <= spi; sti += n, ti++)
                                {
                                    *ti = *sti;
                                }
                                symmetricmatrixvectormultiply(_a, isupper, i + 1, n - 1, t, taui, _t2, n);
                                for (double* st = Tau + i, sp = Tau + n - 2, t2i = _t2; st <= sp; st++, t2i++)
                                {
                                    *st = *t2i;
                                }

                                //
                                // Рассчитываем w := x - 1/2 * tau * (x'*v) * v
                                //
                                v = 0.0d;
                                for (double* st = stj, _taui = Tau + i; st <= spi; st += n, _taui++)
                                {
                                    v += *_taui * *st;
                                }
                                alpha = -(0.5 * taui * v);
                                for (double* st = stj, _taui = Tau + i; st <= spi; st += n, _taui++)
                                {
                                    *_taui += alpha * *st;
                                }

                                //
                                // Применяем преобразование как обновление rank-2:
                                //     A := A - v * w' - w * v'
                                //
                                // 
                                for (double* sti = stj, ti = _t + 1; sti <= spi; sti += n, ti++)
                                {
                                    *ti = *sti;
                                }
                                for (double* st = Tau + i, t2i = _t2 + 1, sp = Tau + n - 2; st <= sp; st++, t2i++)
                                {
                                    *t2i = *st;
                                }
                                symmetricrank2update(_a, isupper, i + 1, n - 1, t, t2, _t3, -1.0d, n);
                                *stj = *E;
                            }
                        *D = *(stj - n);
                        *Tau = taui;
                    }
                    _d[n - 1] = _a[n * n - 1];
                }
        }

        /*************************************************************************
      Генерирует элементарное преобразование отражения

      Подпрограмма генерирует элементарное отражение H порядка N, так, что для
      данного X, имеет место равенство:

          ( X(1) )   ( Beta )
      H * (  ..  ) = (  0   )
          ( X(n) )   (  0   )

      где
                    ( V(1) )
      H = 1 - Tau * (  ..  ) * ( V(1), ..., V(n) )
                    ( V(n) )

      где первая компонента вектора V равна 1.

      Входные параметры:
          X   -   Вектор. Массив с индексами в диапазоне [1..N].
          N   -   Порядок отражения.

      ВЫходные параметры:
          X   -   Кмопоненты от 2 до N, размещенные в векторе V.
                  Первая компонента заменяется с параметром Beta.
          Tau -   Скалярное значение Tau. Если X нулевой вектор, Tau равно 0,
                  иначе 1 <= Tau <= 2.

      Эта подпрограмма является модификацией DLARFG подпрограммы из
      библиотеки LAPACK.

      *************************************************************************/
        static void generatereflection(double* x, int n, out double tau)
        {
            double alpha = 0.0d;
            double xnorm = 0.0d;
            double v = 0.0d;
            double beta = 0.0d;
            double mx = 0.0d;
            double s = 0.0d;

            tau = 0.0d;

            if (n <= 1)
                return;

            //
            // Если требуется, масштабируем (во избежание переполнения/недополнения
            // во время промежуточных вычислений).
            //
            for (double* st = x + 1, sp = x + n; st <= sp; st++)
            {
                mx = Math.Max(Math.Abs(*st), mx);
            }
            s = 1.0d;
            if (mx != 0.0d)
            {
                if (mx <= (minrealnumber / machineepsilon))
                {
                    s = minrealnumber / machineepsilon;
                    v = 1 / s;
                    for (double* st = x + 1, sp = x + n; st <= sp; st++)
                    {
                        *st = v * *st;
                    }
                    mx = mx * v;
                }
                else
                {
                    if (mx >= (maxrealnumber * machineepsilon))
                    {
                        s = maxrealnumber * machineepsilon;
                        v = 1.0d / s;
                        for (double* st = x + 1, sp = x + n; st <= sp; st++)
                        {
                            *st = v * *st;
                        }
                        mx = mx * v;
                    }
                }
            }

            //
            // XNORM = DNRM2( N-1, X, INCX )
            //
            alpha = x[1];
            if (mx != 0.0d)
            {
                for (double* st = x + 2, sp = x + n; st <= sp; st++)
                {
                    xnorm += Sqr(*st / mx);
                }
                xnorm = Math.Sqrt(xnorm) * mx;
            }
            if (xnorm == 0.0d)
            {

                //
                // H  =  I
                //
                x[1] = x[1] * s;
                return;
            }

            //
            // Основной случай
            //
            mx = Math.Max(Math.Abs(alpha), Math.Abs(xnorm));
            beta = -(mx * Math.Sqrt(Sqr(alpha / mx) + Sqr(xnorm / mx)));
            if (alpha < 0.0d)
            {
                beta = -beta;
            }
            tau = (beta - alpha) / beta;
            v = 1.0d / (alpha - beta);
            for (double* st = x + 2, sp = x + n; st <= sp; st++)
            {
                *st = v * *st;
            }
            x[1] = beta;

            //
            // Обратное масштабирование выхода
            //
            x[1] = x[1] * s;
        }

        static void symmetricmatrixvectormultiply(double* a, bool isupper, int i1, int i2,
                                                         double[] X, double alpha, double* y, int N)
        {
            int n = i2 - i1 + 1;

            if (n <= 0)
            {
                return;
            }

            //
            // Примем A = L + D + U, где
            //  L - строго нижний треугольник (нулевая главная диагональ)
            //  D - диагональ
            //  U - строго верхний треугольник (нулевая главная диагональ)
            //
            // A*x = L*x + D*x + U*x
            //
            // Calculate D*x first
            //          
            fixed (double* x = X)
                for (double* st = a + N * i1 + i1, sp = a + N * i2 + i2, _x = x + 1, _y = y + 1; st <= sp; st += N + 1, _x++, _y++)
                {
                    *_y = *st * *_x;
                }

            //
            // Добавляем L*x + U*x
            //
            fixed (double* x = X)
                if (isupper)
                {
                    for (int i = i1; i < i2; i++)
                    {

                        //
                        // Добавляем L*x к результату
                        //
                        double v = x[i - i1 + 1];
                        for (double* st = a + i * N + i + 1, sp = st + n - i - 2 + i1, _y = y + i - i1 + 2; st <= sp; st++, _y++)
                        {
                            *_y = *_y + v * *st;
                        }

                        //
                        // Добавляем U*x к результату
                        //
                        v = 0.0d;
                        for (double* st = a + i * N + i + 1, sp = st + n - i - 2 + i1, _x = x + i - i1 + 2; st <= sp; st++, _x++)
                        {
                            v += *_x * *st;
                        }
                        y[i - i1 + 1] += v;
                    }
                }
                else
                {
                    for (int i = i1 + 1; i <= i2; i++)
                    {

                        //
                        // Добавляем L*x к результату
                        //
                        double v = 0.0d;
                        for (double* _x = x + 1, st = a + i * N + i1, sp = st + i - 1 - i1; st <= sp; st++, _x++)
                        {
                            v += *_x * *st;
                        }
                        y[i - i1 + 1] = y[i - i1 + 1] + v;

                        //
                        // Добавляем U*x к результату
                        //
                        v = x[i - i1 + 1];
                        for (double* _y = y + 1, st = a + i * N + i1, sp = st + i - 1 - i1; st <= sp; st++, _y++)
                        {
                            *_y = *_y + v * *st;
                        }
                    }
                }

            for (double* st = y + 1, sp = y + n; st <= sp; st++)
            {
                *st = alpha * *st;
            }
        }

        static void symmetricrank2update(double* a, bool isupper, int i1, int i2,
                                                double[] X, double[] Y, double* t, double alpha, int N)
        {
            int i = 0;

            fixed (double* x = X, y = Y)
                if (isupper)
                {
                    for (i = i1; i <= i2; i++)
                    {
                        double v = x[i + 1 - i1];
                        for (double* st = y + i + 1 - i1, _t = t + i + 1 - i1, sp = y + i2 - i1 + 1; st <= sp; st++, _t++)
                        {
                            *_t = v * *st;
                        }

                        v = y[i + 1 - i1];
                        for (double* st = x + i + 1 - i1, _t = t + i + 1 - i1, sp = x + i2 - i1 + 1; st <= sp; st++, _t++)
                        {
                            *_t = *_t + v * *st;
                        }

                        for (double* st = t + i + 1 - i1, sp = t + i2 - i1 + 1; st <= sp; st++)
                        {
                            *st = alpha * *st;
                        }

                        for (double* st = a + N * i + i, sp = st + i2 - i, _t = t + i + 1 - i1; st <= sp; _t++, st++)
                        {
                            *st += *_t;
                        }
                    }
                }
                else
                {
                    for (i = i1; i <= i2; i++)
                    {
                        double v = x[i + 1 - i1];

                        for (double* st = y + 1, _t = t + 1, sp = y + i - i1 + 1; st <= sp; st++, _t++)
                        {
                            *_t = v * *st;
                        }

                        v = y[i + 1 - i1];
                        for (double* st = x + 1, _t = t + 1, sp = x + i - i1 + 1; st <= sp; st++, _t++)
                        {
                            *_t += v * *st;
                        }

                        for (double* st = t + 1, sp = t + i - i1 + 1; st <= sp; st++)
                        {
                            *st = alpha * *st;
                        }

                        for (double* st = a + N * i + i1, sp = a + N * i + i, _t = t + 1; st <= sp; st++, _t++)
                        {
                            *st += *_t;
                        }
                    }
                }
        }



        /*************************************************************************
       Распаковка симметричной матрицы Q из трехдиагональной формы.

       Входные параметры:
           A       -   Результат подпрограммы SMatrixTD.
           N       -   Размерность матрицы A.
           IsUpper -   Формат хранения (параметр подпрограммы SMatrixTD).
           Tau     -   результат подпрограммы SMatrixTD.

       Выходные параметры:
           Q       -   Матрица преобразования.
                       Массив с элементами. [0..N-1, 0..N-1].

       *************************************************************************/
        public static void smatrixtdunpackq(double[,] a, int n, bool isupper, double[] tau, ref double[,] q)
        {
            double[] v = new double[0];
            double[] work = new double[0];

            q = new double[n, n];

            if (n == 0)
            {
                return;
            }

            //
            // Инициализация
            //

            v = new double[n + 1];
            work = new double[n];

            fixed (double* _q = q)
                for (double* st = _q, sp = _q + n * n - 1; st <= sp; st += n + 1)
                {
                    *st = 1.0d;
                }

            //
            // Распаковка Q
            //
            fixed (double* w = work, _q = q, _v = v, _a = a)
                if (isupper)
                {
                    for (int j = 0; j < n - 1; j++)
                    {

                        //
                        // Применяем H(i)
                        //
                        for (int i = 0; i < j; i++)
                        {
                            _v[i + 1] = _a[i * n + j + 1];
                        }
                        _v[j + 1] = 1.0d;
                        applyreflectionfromtheleft(_q, tau[j], _v, 0, j, 0, n - 1, w);
                    }
                }
                else
                {
                    for (int j = n - 2; j >= 0; j--)
                    {

                        //
                        // Применяем H(i)
                        //
                        for (int i = 1; i <= n - j - 1; i++)
                        {
                            _v[i] = _a[(i + j) * n + j];
                        }
                        _v[1] = 1.0d;
                        applyreflectionfromtheleft(_q, tau[j], _v, j + 1, n - 1, 0, n - 1, w);
                    }
                }
        }

        /*************************************************************************
        Применение элементарного отражения к прямоугольной матрице размераe MxN.

        Алгоритм предварительно умножает матрицу на преобразование элементарного
        отражения, который получен столбцом V и скаляром Tau (см. описание процедуры
        GenerateReflection). Преобразуется лишь часть матрицы (строки от M1 до M2,
        столбцы от N1 до N2), а не вся матрица. Меняются только элементы данной
        подматрицы.

        Входные параметры:
            C       -   Матрица для преобразования.
            Tau     -   Скаляр определяющий преобразование.
            V       -   Столбец, определяющий преобразование.
                        Массив с индексами в диапазоне [1..M2-M1+1].
            M1, M2  -   Диапазон строк для преобразования.
            N1, N2  -   Диапазон столбцов для преобразования.
            WORK    -   Рабочий массив, индексы которого лежат в диапазоне от N1 до N2.

        Выходные параметры:
            C       -   Результат умножения входной матрицы C с помощью преобразующей
                        матрицы полученной из Tau и V.
                        Если N1>N2 или M1>M2, C не меняется.

        *************************************************************************/
        static void applyreflectionfromtheleft(double* c, double tau, double* v, int m1,
                                                      int m2, int n1, int n2, double* work)
        {
            double t = 0.0d;

            if ((tau == 0.0d || n1 > n2) || m1 > m2)
            {
                return;
            }

            //
            // w := C' * v
            //
            for (double* st = work + n1, sp = work + n2; st <= sp; st++)
            {
                *st = 0.0d;
            }


            for (double* sti = c + (n2 + 1) * m1, spi = c + (n2 + 1) * m2, _v = v + 1; sti <= spi; sti += n2 + 1, _v++)
            {
                for (double* stj = sti + n1, spj = sti + n2, w = work + n1; stj <= spj; stj++, w++)
                {
                    *w += *_v * *stj;
                }
            }

            //
            // C := C - tau * v * w'
            //
            for (double* _v = v + 1, sti = c + (n2 + 1) * m1, spi = c + (n2 + 1) * m2; sti <= spi; sti += n2 + 1, _v++)
            {
                t = *_v * tau;
                for (double* stj = sti + n1, spj = sti + n2, w = work + n1; stj <= spj; stj++, w++)
                {
                    *stj = *stj - t * *w;
                }
            }
        }

        static double Sqr(double p)
        {
            return p * p;
        }

    }  
}

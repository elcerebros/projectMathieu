// zeroin.cpp
/* File converted with FORTRAN CONVERTER utility.
FORTRAN CONVERTER is written by Grigoriev D., 2081/4.*/

#include "FORSYTHE.h"
#include "zeroin.h"

double zeroin(double (*F)(double), double a, double b, double tol) {
    double c, d, e, fa = F(a), fb = F(b), fc, tol1, xm, p, q, r, s;
    while (1) {
        c = a;
        fc = fa;
        e = d = b - a;
        do {
            if (ABS(fc) < ABS(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
            //Проверка сходимости
            tol1 = (double) 2 * EPSILON * ABS(b) + tol / (double) 2;
            xm = (c - b) / (double) 2;
            if (ABS(xm) <= tol1 || fb == 0) return b;
            //Необходима ли бисекция
            if (ABS(e) >= tol1 && ABS(fa) > ABS(fb)) {
                //Возможна ли квадратичная интерполяция
                if (a != c) {
                    //Обратная квадратичная интерполяция
                    q = fa / fc;
                    r = fb / fc;
                    s = fb / fa;
                    p = s * ((double) 2 * xm * q * (q - r) - (b - a) * (r - (double) 1));
                    q = (q - (double) 1) * (r - (double) 1) * (s - (double) 1);
                } else {
                    //Линейная интерполяция
                    s = fb / fa;
                    p = (double) 2 * xm * s;
                    q = (double) 1 - s;
                }
                //Выбрать знаки
                if (p > 0) q = -q;
                else p = -p;//p = ABS(p)

                //Приемлема ли интерполяция
                if (((double) 2 * p) >= ((double) 3 * xm * q - ABS(tol1 * q)) || p >= ABS(e * q / (double) 2))
                    e = d = xm;//Бисекция
                else {
                    e = d;
                    d = p / q;
                }
            } else e = d = xm;//Бисекция
            //Завершить шаг
            a = b;
            fa = fb;
            if (ABS(d) > tol1) b += d;
            else b += SIGN2(tol1, xm);
            fb = F(b);
        } while (fb * SIGN(fc) <= 0);
    }
}


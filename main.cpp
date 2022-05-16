#include <iostream>
#include <cmath>
#include "quanc8.h"
#include "zeroin.h"
#include "rkf45.h"

using namespace std;

// Инициализация переменных и констант
double E;
double A;
double delta;
const int N = 20;
const float B = 0.0;

// Изначальное дифференциальное уравнение, приведёное к системе уравнений первого порядка
void func(double t, const double *y, double* dy) {
    dy[0] = y[1];
    dy[1] = -y[0] * (delta - E * cos(2 * t));
}

// Подынтегральная функция при вычислении константы E
double funcComputingE(double x) {
    return sin(x) / (pow(x, 2) + 1);
}

// Функция, необходимая для вычисления константы A
double funcComputingA(double x) {
    return x - pow(1.4, x);
}

double ce2n(double x) {
    // 1 * x = 0.25 * y
    // y = x / 0.25

    // 1 0.25
    // 0.25 0.0625
    // 0.625

    double result = 1 * cos(0 * x) + 0.25 * cos(2 * x) + 0.0625 * cos(4 * x) + 0.015625 * cos(6 * x);
    cout << fixed << "x = " << x << "\ty = " << result << endl;

    return result;
}

double quanc8Computing() {
    double abserr = 0.0;
    double relerr = 0.0;
    double errest, flag;
    int nofun;
    double left = 0.0;
    double right = 1.0;
    double resultQuanc;

    quanc8(funcComputingE, left, right, abserr, relerr, &resultQuanc, &errest,
           &nofun, &flag);

    cout << "\nQUANC8 method:" << endl;
    cout << "\tValue of the integrand (in E): " << resultQuanc << endl;
    cout << "\tError: " << errest << endl;

    return resultQuanc;
}

double zeroinComputing() {
    double a = 0.0;
    double b = 4.0;
    double eps = 1e-8;

    double result = zeroin(funcComputingA, a, b, eps);

    cout << "\nZEROIN method:" << endl;
    cout << "\tThe root of the equation: " << result << endl;
    cout << "\tError: " << eps << endl;

    return result;
}

void rkf45(double* u, double* du, double u0, double du0, double h) {
    int iflag = 1;
    int iwork[30];
    double t = 0.0;
    double tout = 0.0;
    double work[15];
    double RE = 1e-10;
    double AE = 1e-10;
    double U0[] = {u0, du0};

    cout << "\tRKF45 method:" << endl;
    cout.precision(10);

    for (int i = 0; i < N; ++i) {
        RKF45(func, 2, U0, &t, &tout, &RE, &AE, &iflag,
              work, iwork);

        cout << fixed << "\ttout = " << tout << "\t" << left << "u = " << U0[0] << "\t" << "du = " << U0[1]
             << "\tiflag = " << iflag << endl;

        u[i] = U0[0];
        du[i] = U0[1];
        tout += h;
    }
}

int main() {
    double quancResult = quanc8Computing();
    double zeroinResult = zeroinComputing();

    // Вычисление констант, необходимых для решения уравления Матье
    E = 1.553791 * quancResult;
    A = 0.5300355 * zeroinResult;
    cout << "\nCOMPUTING OF CONSTANTS:" << endl;
    cout << "\t1) E = " << E << endl;
    cout << "\t2) A = " << A << endl;

    // Решение уравнения Матье
    double rkfU0[N], rkfU1[N];
    float U0[] = {static_cast<float>(A), B};
    float h = 0.5;
    cout << "\nSOLUTION:" << endl;
    rkf45(rkfU0, rkfU1, U0[0], U0[1], h);

    cout << "\nSOLUTION REAL:" << endl;
    double x = 0.0;
    while (x <= 10.0) {
        ce2n(x);
        x += 0.5;
    }

    return 0;
}
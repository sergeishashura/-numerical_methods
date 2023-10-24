#include <iomanip>
#include "headerFunctions.h"

// Функция для решения системы уравнений Ax = b методом Гаусса
// A - матрица коэффициентов размера 3x3
// b - вектор правой части размера 3
// x - вектор решения размера 3

void pause() {
    cout << "\npress enter to continue... \n";
    cin.clear();
    cin.sync();
    cin.get();
}

vector<vector<double>> initA() {
    vector<vector<double>> A = {{2.31, 31.49,  1.52},
                                {4.21, 22.42, 3.85},
                                {3.49, 4.85,  28.72}};
    return A;
}

vector<double> initB() {
    vector<double> A = {40.95, 30.24, 42.81};
    return A;
}

vector<double> countResidualVector(vector<vector<double>> A, vector<double> b, vector<double> x) {
    std::vector<double> F(3, 0.0);
    int n = F.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            F[i] += A[i][j] * x[j];
        }
    }
//    vector<double> tmp(3, 0.0);
    for (int i = 0; i < n; ++i) {
        F[i] = F[i] - b[i];
    }
    return F;
}

vector<double> solveAuxiliary(vector<vector<double>> A, vector<double> x) {
    int n = 3;
    vector<double> x1(n);
    vector<double> x11(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            x1[i] += A[i][j] * x[j];
        }
    }
    x11 = solveGauss(A, x1);
    return x11;
}

double countRelativeError(vector<double> x, vector<double> x1) {
    double norm1 = 0;
    double norm2 = 0;
    double error = 0;

    for (int i = 0; i < x.size(); i++) {
        norm1 += (x[i] - x1[i]) * (x[i] - x1[i]);
        norm2 += x1[i] * x1[i];
    }
    norm1 = sqrt(norm1);
    norm2 = sqrt(norm2);
    error = norm1 / norm2;

    return error;
}


vector<double> solveGauss(vector<vector<double>> A, vector<double> b) {
    int n = 3;
    vector<double> x(n);
    // прямой ход - приводим мтрицу к нужному виду
    for (int k = 0; k < n - 1; k++) {
        for (int i = k + 1; i < n; i++) {

            double m = A[i][k] / A[k][k]; // коэффициент, на который умножаем k-ую строку
            for (int j = k; j < n; j++) {
                A[i][j] = A[i][j] - m *
                                    A[k][j]; // вычитаем из элемента i-ой строки произведение коэффициента и элемента k-ой строки

            }
            b[i] = b[i] - m * b[k]; // аналагично для вектора правой части
        }
        output(A, b, x);
        pause();
        for (signed i = k; i < n; i++) {
            double d = A[i][i]; // диагональный элемент
            for (int j = k; j < n; j++) {
                A[i][j] = A[i][j] / d; // делим элемент на диагональный
            }
            b[i] = b[i] / d;
        }
    }
    output(A, b, x);
    pause();


// обратный ход - находим решения по формуле x[i] = b[i] - sum(A[i][j] * x[j]) где j > i
    for (int i = n - 1; i >= 0; i--) {
        double s = 0; // сумма произведений коэффициентов и уже найденных решений
        for (int j = i + 1; j < n; j++) {
            s += A[i][j] * x[j]; // прибавляем к сумме произведение коэффициента и решения
        }
        x[i] = b[i] - s;
    }
    cout << "\n" << "After counting: \n";
    output(A, b, x);
    return x;
}

double countShit(vector<vector<double>> A, vector<double> b, vector<double> x) {
    vector<double> r = countResidualVector(A, b, x);
    double norm = 0;
    for (int i = 0; i < r.size(); i++) {
        norm += r[i] * r[i];
    }
    norm = sqrt(norm);
    return norm;
}

void output(vector<vector<double>> A, vector<double> b, vector<double> x) {

    cout << "Matrix A: \n";
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].size(); ++j) {
            cout << A[i][j] << setw(20);
        }
        cout << "\t |" << b[i] << "\n";
    }

    cout << "\n X:";
    for (int i = 0; i < x.size(); ++i) {
        cout << "x" << i + 1 << " =" << x[i] << "\t";
    }
}
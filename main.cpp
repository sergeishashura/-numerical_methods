#include "headerFunctions.h"

int main() {
    vector<vector<double>> A = initA();
    vector<double> b = initB();
    vector<double> x(3,0.0);
    output(A,b,x);
    x = solveGauss(A, b);
    vector<double> F = countResidualVector(initA(), initB(), x);
    vector<double> x11 = solveAuxiliary(A, x);
    double shit = countRelativeError(x,x11);
    cout << "\nVector F:\t";
    for(double i  : F) {
        cout << i << "\t";
    }
    cout << "\nShit: \t" << shit;
}
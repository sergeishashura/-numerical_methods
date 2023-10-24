#ifndef LABA_1_HEADERFUNCTIONS_H
#define LABA_1_HEADERFUNCTIONS_H

#endif //LABA_1_HEADERFUNCTIONS_H

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
vector<vector<double>> initA();
vector<double> initB();
double countShit(vector<vector<double>> A, vector<double> b, vector<double> x);
vector<double> solveAuxiliary(vector<vector<double>> A, vector<double> x);
double countRelativeError(vector<double> x, vector<double> x1);
vector<double>  countResidualVector (vector<vector<double>> A, vector<double> b, vector<double> x);
vector<double> solveGauss(vector<vector<double>> A, vector<double> b);
void output(vector<vector<double>> A, vector<double> b, vector<double> x);
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

double MultiplyParametersOfVectors(vector<double> firstVector, vector<double> secondVector) // Σ(fV[i] * sV[i])
{
    double resultOfMult = 0;
    for (int i = 0; i < firstVector.size(); i++)
        resultOfMult += firstVector[i] * secondVector[i];
    return resultOfMult;
}

double SumVectorElements(vector<double> vectorOfElements)
{
    double sumOfVectorElements = 0;
    for (int i = 0; i < vectorOfElements.size(); i++)
        sumOfVectorElements += vectorOfElements[i];
    return sumOfVectorElements;
}

vector<double> RaiseElementsOfVectorToPow(vector<double> vectorOfElements, int degree)
{
    vector<double> vectorOfResults(vectorOfElements.size());
    for (int i = 0; i < vectorOfElements.size(); i++)
        vectorOfResults[i] = pow(vectorOfElements[i], degree);
    return vectorOfResults;
}

vector<double> SolveLSMApproximation(vector<double> experimentalValuesOfX, vector<double> experimentalValuesOfY)
{
    int numberOfExperimentalPoints = experimentalValuesOfX.size();
    int polynomialDegree = 2; // Quadratic polynomial

    vector<vector<double>> coefficientMatrix(polynomialDegree + 1, vector<double>(polynomialDegree + 1, 0.0));
    vector<double> vectorOfMultiplyParameters(polynomialDegree + 1, 0.0);

    for (int i = 0; i <= polynomialDegree; i++)
    {
        for (int j = 0; j <= polynomialDegree; j++)
        {
            if (i != 0 || j != 0)
            {
                int k = i + j;
                coefficientMatrix[i][j] = SumVectorElements(RaiseElementsOfVectorToPow(experimentalValuesOfX, k));
            }
            else
                coefficientMatrix[i][j] = numberOfExperimentalPoints;
        }
        vectorOfMultiplyParameters[i] = MultiplyParametersOfVectors(experimentalValuesOfY, RaiseElementsOfVectorToPow(experimentalValuesOfX, i));
    }

    vector<double> equationCoefficients = solveGauss(coefficientMatrix, vectorOfMultiplyParameters);
    return equationCoefficients;
}

double CalculateStandardDeviation(vector<double> experimentalValuesOfX, vector<double> experimentalValuesOfY, vector<double> equationCoefficients)
{
    double standardDeviation = 0;
    int numberOfExperimentalPoints = experimentalValuesOfX.size();
    int polynomialDegree = equationCoefficients.size() - 1;

    for (int i = 0; i < numberOfExperimentalPoints; i++)
    {
        double sumOfVectorElements = 0;
        for (int j = 0; j <= polynomialDegree; j++)
            sumOfVectorElements += equationCoefficients[j] * pow(experimentalValuesOfX[i], j);
        standardDeviation += pow(experimentalValuesOfY[i] - sumOfVectorElements, 2);
    }

    standardDeviation /= (numberOfExperimentalPoints - polynomialDegree - 1);
    return standardDeviation;
}

int main()
{
    vector<double> experimentalValuesOfX = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    vector<double> experimentalValuesOfY = {0.957, 0.969, 0.976, 0.978, 0.975, 0.968, 0.954, 0.939, 0.918, 0.894};

    vector<double> equationCoefficients = SolveLSMApproximation(experimentalValuesOfX, experimentalValuesOfY);
    cout << "Quadratic polynomial equation: v = a + bD + cD^2\n";
    cout << "Found equation coefficients:\n";
    cout << "a = " << equationCoefficients[0] << ", b = " << equationCoefficients[1] << ", c = " << equationCoefficients[2] << endl;

    double standardDeviation = CalculateStandardDeviation(experimentalValuesOfX, experimentalValuesOfY, equationCoefficients);
    cout << "Standard deviation of the equation: " << standardDeviation << endl;

    return 0;
}
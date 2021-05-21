/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Sooyeon Han
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"



Matrix   linearFit(Matrix _x, Matrix _y) {

    int mx = _x.rows;
    int my = _y.rows;

    double a1 = 0;
    double a0 = 0;

    int n = mx;

    double Sxx = 0;
    double Sx = 0;
    double Sxy = 0;
    double Sy = 0;

    for (int i = 0; i < n; i++)
    {
        Sxx += pow(_x.at[i][0], 2);
        Sx += _x.at[i][0];
        Sxy += _x.at[i][0] * _y.at[i][0];
        Sy += _y.at[i][0];
    }

    a1 = (n * Sxy - Sx * Sy) / (n * Sxx - Sx * Sx);
    a0 = (Sxx * Sy - Sxy * Sx) / (n * Sxx - Sx * Sx);

    double z[2] = { a1, a0 };

    for (int i = 0; i < 2; i++)
    {
        printf("%d\n", z[i]);
    }

    double z_array[] = { a1, a0 };
    return arr2Mat(z_array, 2, 1);
}



// Create a matrix from 1D-array

Matrix   arr2Mat(double* _1Darray, int _rows, int _cols)

{
    Matrix Output = createMat(_rows, _cols);
    for (int i = 0; i < _rows; i++)
        for (int j = 0; j < _cols; j++)
            Output.at[i][j] = _1Darray[i * _cols + j];
    return Output;

}



double linearInterp(Matrix _x, Matrix _y, double xq)

{
    double num = _x.rows;
    double yq = 0;

    int i = 0;
    for (i = 0; i < num; i++)
    {
        // i를 찾기 위해 if문을 써 준다
        if (xq < _x.at[i][0])
        {
            break;
        }
    }

    yq = _y.at[i - 1][0] * (xq - _x.at[i][0]) / (_x.at[i - 1][0] - _x.at[i][0]) + _y.at[i][0] * (xq - _x.at[i - 1][0]) / (_x.at[i][0] - _x.at[i - 1][0]);
    return yq;

}


Matrix	gradient(Matrix _x, Matrix _y)
{
    double h = 0;
    Matrix f = createMat(_x.rows, 1);

    if (_x.rows != _y.rows)
    {
        printf("error!!\n size is not same.");
        return f;
    }

    if (_x.rows == 2)
    {
        h = _x.at[1][0] - _x.at[0][0];
        f.at[0][0] = (_y.at[1][0] - _y.at[0][0]) / h;
        return f;
    }


    for (int i = 0; i < _x.rows; i++)
    {
        if (i == 0)
        {
            h = _x.at[i + 1][0] - _x.at[i][0];
            f.at[i][0] = (-3 * _y.at[i][0] + 4 * _y.at[i + 1][0] - _y.at[i + 2][0]) / (2 * h);

        }
        else if (i == _x.rows - 1)
        {
            h = _x.at[i][0] - _x.at[i - 1][0];
            f.at[i][0] = (_y.at[i - 2][0] - 4 * _y.at[i - 1][0] + 3 * _y.at[i][0]) / (2 * h);
        }

        else
        {
            h = _x.at[i + 1][0] - _x.at[i - 1][0];
            f.at[i][0] = (_y.at[i + 1][0] - _y.at[i - 1][0]) / h;
        }
    }
    return f;
}


void gradient1D(double _x[], double _y[], double _dydx[], int m)
{
    double h = 0;

    for (int i = 0; i < m; i++)
    {
        if (i == 0)
        {
            h = _x[i + 1] - _x[i];
            _dydx[i] = (-3 * _y[i] + 4 * _y[i + 1] - _y[i + 2]) / (2 * h);
            //printf("hh%f\n", f.at[i][0]);

        }
        else if (i == m - 1)
        {
            h = _x[i] - _x[i - 1];
            _dydx[i] = (_y[i - 2] - 4 * _y[i - 1] + 3 * _y[i]) / (2 * h);
        }

        else
        {
            h = _x[i + 1] - _x[i - 1];
            _dydx[i] = (_y[i + 1] - _y[i - 1]) / h ;
        }
    }
}

// Define a function that defines the target equation.
double myFunc(const double x) {
    return  x * x * x;
}

// Return the dy/dx results for the target equation. (truncation error: O(h^2))
// Move this function to myNM.h and myNM.cpp
Matrix	gradientFunc(double func(const double x), Matrix _xin) 
{

    int n = _xin.rows;
    Matrix y = createMat(n, 1);
    for (int i = 0; i < _xin.rows; i++)
    {
        y.at[i][0] = func(_xin.at[i][0]);
    }
    y = gradient(_xin, y);

    return y;
}

// Define a function that defines the target equation.
double dmyFunc(const double x) 
{
    return  3 * x * x;
}

double newtonRaphsonFunc(double func(const double x), double dfunc(const double x), double _x0, double tol)
{
    double x = _x0;
    int N_max = 20;
    double _h = 0;
    int N = 0;

    do {
        _h = -func(x) / dfunc(x);
        x = x + _h;

        printf("Iteration:%d \t", N);
        printf("X(n): %f \t", x);
        printf("Tolerance: %.10f\n", fabs(func(x)));
        N++;
    } while (N < N_max && fabs(func(x)) >= tol);

    return x;
}


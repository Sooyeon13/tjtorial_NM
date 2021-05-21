/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Sooyeon Han
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H

#include "myMatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern Matrix linearFit(Matrix _x, Matrix _y);

extern Matrix arr2Mat(double* _1Darray, int _rows, int _cols);

extern double linearInterp(Matrix _x, Matrix _y, double xq);

extern Matrix gradient(Matrix _x, Matrix _y);

extern void gradient1D(double x[], double y[], double dydx[], int m);

extern Matrix gradientFunc(double func(const double x), Matrix xin);

extern double myFunc(const double x);

extern Matrix gradientFunc(double func(const double x), Matrix _xin);

extern double dmyFunc(const double x);

extern double newtonRaphsonFunc(double func(const double x), double dfunc(const double x), double x0, double tol);
#endif
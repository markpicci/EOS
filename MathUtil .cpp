#include "MathUtil.h"
#include <format>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cmath>

//==============================================================================================
//simple matrix dump
//
void matrixDump(std::vector<std::vector<double>>& m){
    std::filesystem::path directorypath;
    directorypath = std::filesystem::current_path();
    std::string fileName = directorypath.string();
    fileName += "/matrix.txt";
    std::ofstream outputFileStream(fileName, std::ios::out);

    for(const std::vector<double>& r : m){
        for(const double& d : r){
            if (std::isnan(d)) outputFileStream << std::format("{:<8s} ","NAN");
            else outputFileStream << std::format("{:<8.4f} ",d);
        }
        outputFileStream << "\n";
    }
    outputFileStream.close();
}

bool checkSymmetric(std::vector<std::vector<double>>& m){
    size_t row = m.size();
    size_t col = m[0].size();

    if (row!=col) return false;
    for(size_t rowi=0; rowi<row; rowi ++)
        for(size_t coli=0; coli<rowi; coli ++)
            if (m[rowi][coli]!=m[coli][rowi]) return false;
    return true;
}

// Kahan summation algorithm
// Function to implement the Kahan
// summation algorithm
double kahanSum(std::vector<double> &fa)
{
    double sum = 0.0;

    // Variable to store the error
    double c = 0.0;

    // Loop to iterate over the array
    for(double f : fa)
    {
        double y = f - c;
        double t = sum + y;

        // Algebraically, c is always 0
        // when t is replaced by its
        // value from the above expression.
        // But, when there is a loss,
        // the higher-order y is cancelled
        // out by subtracting y from c and
        // all that remains is the
        // lower-order error in c
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

//==============================================================================================
//simple bubble sort
//
template <typename T>
void bubbleSort(T arr[], int n)
{
    int i, j;
    for (i = 0; i < n - 1; i++)

        // Last i elements are already
        // in place
        for (j = 0; j < n - i - 1; j++)
            if (arr[j] > arr[j + 1])
                std::swap(arr[j], arr[j + 1]);
}

//==============================================================================================
//function for temperature conversion.
//
double Tcon(char in, double val, char out){
double conv;

    conv =0.0;
    if(in=='C' and out=='F') conv = val*9.0/5.0+32.0;
    if(in=='F' and out=='C') conv = (val-32.0)*5.0/9.0;

    if(in=='C' and out=='K') conv = val+273.15;
    if(in=='R' and out=='K') conv = val*5.0/9.0;
    if(in=='F' and out=='K') conv = 5.0/9.0*(val+459.67);

    if(in=='K' and out=='C') conv = val-273.15;
    if(in=='K' and out=='R') conv = val*9.0/5.0;         //Temperature 1.8 °R = 1.0 K (exact)
    if(in=='K' and out=='F') conv = val*9.0/5.0-459.67;
    return conv;
}

//==============================================================================================
//estimate the root of cubic expression using Cardano's method.
//
double Xroot(double a, double x)
{
    double i = 1;
    if (a < 0) i = -1;
    return (i * exp(log(a*i)/x));
}

int cubicRoot(double a1, double b1, double c1, double d1, double root[], double& img)  // solve cubic equation according to cardano
{
    double a, b, c;
    double p, q, u, v, d;
    double r, alpha;
    double x1real, x2real, x2imag, x3real, x3imag;
    int res;
    res = 0;
    if (a1 != 0)
    {
        a = b1 / a1;
        b = c1 / a1;
        c = d1 / a1;

        p = -(a * a / 3.0) + b;
        q = (2.0 / 27.0 * a * a * a) - (a * b / 3.0) + c;
        d = q * q / 4.0 + p * p * p / 27.0;
        // 3 cases D > 0, D == 0 and D < 0
        if (d > epsD)
        {
            u = Xroot(-q / 2.0 + sqrt(d), 3.0);
            v = Xroot(-q / 2.0 - sqrt(d), 3.0);
            x1real = u + v - a / 3.0;
            x2real = -(u + v) / 2.0 - a / 3.0;
            x2imag = sqrt(3.0) / 2.0 * (u - v);
            x3real = x2real;
            x3imag = -x2imag;
            root[0]=x1real;root[1]=x2real;root[2]=x3real;
            img = x2imag;
            res = 1;
        }
        if (abs(d) <= epsD)
        {
            u = Xroot(-q / 2.0, 3.0);
            v = Xroot(-q / 2.0, 3.0);
            x1real = u + v - a / 3.0;
            x2real = -(u + v) / 2.0 - a / 3.0;
            root[0]=x1real;root[1]=x2real;root[2]=0.0;
            img = 0.0;
            res = 2;
            if(root[0] > root[1]) std::swap(root[0], root[1]);
        }
        else
        {
            r = sqrt(-p * p * p / 27.0);
            alpha = atan(sqrt(-d) / (-q) * 2.0);
            if (q > 0)                         // if q > 0 the angle becomes  PI + alpha
                alpha = pi + alpha;

            x1real = Xroot(r, 3.0) * (cos((6.0 * pi - alpha) / 3.0) + cos(alpha / 3.0)) - a / 3.0;
            x2real = Xroot(r, 3.0) * (cos((2.0 * pi + alpha) / 3.0) + cos((4.0 * pi - alpha) / 3.0)) - a / 3.0;
            x3real = Xroot(r, 3.0) * (cos((4.0 * pi + alpha) / 3.0) + cos((2.0 * pi - alpha) / 3.0)) - a / 3.0;
            root[0]=x1real;root[1]=x2real;root[2]=x3real;
            img = 0.0;
            res = 3;
            bubbleSort<double>(root,3);
        }
    }

    return res;
 }

//==============================================================================================
//estimate the root using Chamber's method (quadratic solution).
//
int Zero(double& sol, double goal, double fun(double), double X1, double X2) {
long I;
double X3, F, F1, F2, F3;
double TOL = 0.5e-9;

    sol = 0.0;
    F1=fun(X1)-goal;
    F2=fun(X2)-goal;
    if (F1*F2>=0) return 0;

    //------------------
    // BEGIN ITERATING
    //------------------
    for (I=1; I<=50; I++){
    // Use False Position to get point 3.*/
        X3 = X1-F1*(X2-X1)/(F2-F1);
        F3=fun(X3)-goal;
    // Use points 1, 2, and 3 to estimate the root using Chamber's method (quadratic solution).
        sol = X1*F2*F3/((F1-F2)*(F1-F3)) + X2*F1*F3/((F2-F1)*(F2-F3)) + X3*F1*F2/((F3-F1)*(F3-F2));
        if ((sol-X1)*(sol-X2)>=0) sol=(X1+X2)/2;
        F=fun(sol)-goal;
        if (fabs(F)<=TOL) return I;
    // Discard quadratic solution if false position root is closer.
        if (fabs(F3) < fabs(F) && F*F3 > 0){
            if (F3*F1>0) {X1=X3; F1=F3;}
            else {X2=X3; F2=F3;}}
        else{
    // Swap in new value from quadratic solution
            if (F*F3 < 0) {X1=sol; F1=F; X2=X3; F2=F3;}
            else if (F3*F1 > 0) {X1=sol; F1=F;}
            else {X2=sol; F2=F;}}
    }
    sol = 0.0;
    return I;
}

//==============================================================================================
//estimate the root using Bisection Method.
//
int ZeroBisec(double& sol, double goal, double fun(double), double X1, double X2) {
long itr, maxitr=50;
double F, F1, F2;
double TOL = 0.5e-9;

    sol = 0.0;
    itr = 0;
    F1=fun(X1)-goal;
    F2=fun(X2)-goal;
    if (F1*F2>=0) return itr;

    //------------------
    // BEGIN ITERATING
    //------------------
    for(itr=1;itr<=maxitr;itr++){
        /* Bisecting Interval */
        sol = (X1 + X2)/2.0;
        F = fun(sol)-goal;
        if(fabs(F)<TOL) return itr;
        if(F1*F < 0){
            X2 = sol;
            F2 = F;
        } else {
            X1 = sol;
            F1 = F;
        }
	 }

	 return itr;
}

//==============================================================================================
//estimate the root using Newton Method.
//
int ZeroNew(double& sol, double goal, double fun(double),double funDer(double)) {
int itr, maxitr=50;
double h, TOL = 0.5e-9;

    for(itr=1;itr<=maxitr;itr++){
        h=(fun(sol)-goal)/funDer(sol);
        sol-=h;
        if(fabs(h)<TOL) return itr;
    }
    return itr;
}

//==============================================================================================
// Binary Search [n, n+1)
// CHECK x<xd[low] and x>xd[high] shall be done before enter the sub
// high SAHLL be a valid entry array sized (low, high+1)
template <class T> inline size_t binarySearch(const std::vector<T>& xd, T x, size_t low, size_t high)
{
    size_t left = low;
    size_t half = high-low+1;

    while (half > 1) {
        half >>= 1;
        left = (x < xd[left + half] ? left : left+half);
    }
    return left;
}
//==============================================================================================
// Linear Search [n, n+1)
// CHECK x<xd[low] and x>xd[high] shall be done before enter the sub
// high SAHLL be a valid entry array sized (low, high+1)
template <class T> inline size_t linearSearch(const std::vector<T>& xd, T x, size_t low, size_t high)
{
    for(size_t i = low; i< high; i++)
        if (x>=xd[i] && x<xd[i+1]) return i;
    return high-1;
}

//==============================================================================================
//Linear Interpolation
//
int interLin(double x, double& y, const std::vector<double>& xv, const std::vector<double>& yv, int extrapolate = 1) {
int i, nEle;

//Initial check
    if(xv.size()!=yv.size()) return 1; //ERROR not same size
    nEle = xv.size();
    y =0.0;
    if(nEle<2) return 1; //ERROR we need at least two point

    switch(extrapolate){
        case 0: // NO Extrapolation
            if(x<xv[0])    return 2; //underflow
            if(x>xv[nEle-1])  return 3; //overflow
            break;

        case 1: // Extrapolation with value
            if(x<xv[0])    {y=yv[0];  return 4;}  //underflow
            if(x>xv[nEle-1])  {y=yv[nEle-1];return 5;} //overflow
            break;

        case 2: break;

        default:
            if(x == xv[0]) {y=yv[0]; return 0;}
            if(x == xv[nEle-1]) {y=yv[nEle-1]; return 0;}
    }

    i=binarySearch<double>(xv, x, 0, nEle-1);
    if(i > nEle-2) i = nEle-2;

//Calculate INTERPOLATION in interval i-1 / i
    double m    = yv[i+1]-yv[i];
    double den  = xv[i+1]-xv[i];
    if(fabs(den)< epsD ){
            y = (yv[i+1]+yv[i]) / 2.0; //bad situation
        } else {
            m /= den;
            y =  m* (x-xv[i]) + yv[i];
        }

    return 0;
}

//****************************************************************************80
//
//  Purpose:
//    LAGRANGE_VALUE_1D evaluates the Lagrange interpolant.
//
//  Discussion:
//    The Lagrange interpolant L(ND,XD,YD)(X) is the unique polynomial of
//    degree ND-1 which interpolates the points (XD(I),YD(I)) for I = 1
//    to ND.
//
//    The Lagrange interpolant can be constructed from the Lagrange basis
//    polynomials.  Given ND distinct abscissas, XD(1:ND), the I-th Lagrange
//    basis polynomial LB(ND,XD,I)(X) is defined as the polynomial of degree
//    ND - 1 which is 1 at  XD(I) and 0 at the ND - 1 other abscissas.
//
//    Given data values YD at each of the abscissas, the value of the
//    Lagrange interpolant may be written as
//
//      L(ND,XD,YD)(X) = sum ( 1 <= I <= ND ) LB(ND,XD,I)(X) * YD(I)
//
//  Parameters:
//
//    Input, double XD[ND], the data points.
//    Input, double YD[ND], the data values.
//
//    Output, double LAGRANGE_VALUE_1D, the interpolated values.
//
int lagrange_value_1d (double x, double& y, const std::vector<double>& xd, const std::vector<double>& yd, int extrapolate){
int n, i, j, index, order;
std::vector<double> lb(4); //order  = 4

//Initial check
    order = 4;
    //lb.resize(order)
    if(xd.size()!=yd.size()) return 1; //ERROR not same size
    n = xd.size();
    if (n<=1) return 1; //ERROR we need at least two point
    if (n<=4) order = n; // Order is order of polinomial order 3 require 4 point
    if (n>4)  order = 4; // Order is order of polinomial order 3 require 4 point

    switch(extrapolate){
        case 0: // NO Extrapolation
            if(x<xd[0])    return 2; //underflow
            if(x>xd[n-1])  return 3; //overflow
            break;

        case 1: // Extrapolation with value
            if(x<xd[0])    {y=yd[0];  return 4;}  //underflow
            if(x>xd[n-1])  {y=yd[n-1];return 5;} //overflow
            break;

        case 2: break;

        default:
            if(x == xd[0]) {y=yd[0]; return 0;}
            if(x == xd[n-1]) {y=yd[n-1]; return 0;}
    }

    index = 0;
    // If order is less then 2 linear or parabolic we have 2 or three point no need to look for
    if (n>4){
        index=binarySearch<double>(xd, x, 0, n-1);
        if (index>=n-order) index = n-order;
        if (index>=1 && index<n-order) index--; // keep interpolation on middle span
    }

    //Lagrange Base
    for ( j = 0; j < order; j++ ) lb[j] = 1.0;
    for ( i = 0; i < order; i++ )
        for ( j = 0; j < order; j++ )
            if ( j != i )
            lb[i] *= ( x- xd[j+index] ) / ( xd[i+index] - xd[j+index] );

    y = 0.0;
    for (int j = 0; j < order; j++)
        y += lb[j] * yd[j+index];

    return 0;
}

int nofDecimals(const double val)
{
    int result = 0;
    double epsilon = epsD;
    double exponent = 1.0;

    while(fabs(fmod(val * exponent, 1.0)) > epsilon){
        ++result;
        epsilon  *= 10;
        exponent *= 10;
    }

    return result;
}

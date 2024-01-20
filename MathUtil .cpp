#include "MathUtil.h"
#include <iostream>

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
//estimate the root using Chamber's method (quadratic solution).
//
int ZeroBisec(double& sol, double goal, double fun(double), double X1, double X2) {
long I, itr, maxitr=50;
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
/*
// Dummy for test
double fun(double x){
    //double ret = x*x*x-2.9*x*x+2.8*x-0.9;
    double ret = pow(2,x);
    return ret ;
}
// Dummy for test derivate
double funDer(double x){
    //double ret = x*x*x-2.9*x*x+2.8*x-0.9;
    double ret = log(2.0)*pow(2,x);
    return ret ;
}

main(){
double root[3], img, ris;
int nRoot;
double a, b, c;
double x, x2, x3;

    //a = -4.0, b = 4.0, c = 6.0;  //1
    //a = 4.0, b = -34.0, c = 6.0; //3
    a = -2.9, b = 2.8, c = -0.9;   //2

    nRoot = cubicRoot(1, a, b, c, root, img);
    std::cout << "nroots " << nRoot << "     - EPS "<< epsD << "\n";
    std::cout << root[0] << " "<< root[1] << " "<< root[2] << "\n";
    std::cout << "Img "<< img << "\n\n";

    x = root[0];
    x2 = x*x;
    x3 = x*x*x;
    std::cout << " Zero 1 " << x3+a*x2+b*x+c << "\n";
    if (nRoot==2){
        x = root[1];
        x2 = x*x;
        x3 = x*x*x;
        std::cout << " Zero 2 " << x3+a*x2+b*x+c << "\n";
    }
    if (nRoot==3){
        x = root[1];
        x2 = x*x;
        x3 = x*x*x;
        std::cout << " Zero 2 " << x3+a*x2+b*x+c << "\n";
        x = root[2];
        x2 = x*x;
        x3 = x*x*x;
        std::cout << " Zero 3 " << x3+a*x2+b*x+c << "\n";
    }

    std::cout << "\n";
    int nIter = Zero(x, 5.0, fun, 0.1, 5.0);
    std::cout << " Chamber's method iteration " << nIter << "\n";
    std::cout << " Chamber's method X " << x << "\n";
    std::cout << " Chamber's method Y " << fun(x) << "\n";

    std::cout << "\n";
    nIter = ZeroNew(x, 5.0, fun, funDer);
    std::cout << " Newton's method iteration " << nIter << "\n";
    std::cout << " Newton's method X " << x << "\n";
    std::cout << " Newton's method Y " << fun(x) << "\n";

    std::cout << "\n";
    nIter = ZeroBisec(x, 5.0, fun, 0.1, 5.0);
    std::cout << " Bisec's method iteration " << nIter << "\n";
    std::cout << " Bisec's method X " << x << "\n";
    std::cout << " Bisec's method Y " << fun(x) << "\n";

}*/

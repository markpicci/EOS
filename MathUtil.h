#ifndef _MATH_UTIL_H
#define _MATH_UTIL_H

#define _USE_MATH_DEFINES
#include <cmath>

const double pi         = M_PI;         // pi
const double pi_2       = M_PI_2;       // pi/2
const double pi_4       = M_PI_4;       // pi/4
const double p1_pi      = M_1_PI;       // 1/pi
const double p2_pi      = M_2_PI;       // 2/pi
const double p2_sqrtpi  = M_2_SQRTPI;   // 2/sqrt(pi)
const double sqrt2      = M_SQRT2;      // sqrt(2)
const double sqrt1_2    = M_SQRT1_2;    // 1/sqrt(2)
const double e          = M_E;          // e
const double log2e      = M_LOG2E;      // log_2(e)
const double log10e     = M_LOG10E;     // log_10(e)
const double ln2        = M_LN2;        // log_e(2)
const double ln10       = M_LN10;       // log_e(10)
const double epsD       = std::numeric_limits<double>::epsilon();
const double inf        = std::numeric_limits<double>::infinity();
//const double nan        = std::numeric_limits<double>::quiet_NaN()
const double g          = 9.806650;      //m/s
const double gUS        = 32.17405;      //ft/s

const long double piL = 3.141592653589793238462643383279502884L; /* pi */

const double incTOm  = 0.0254; //Length 1.0 in = 0.0254 m  (exact)
const double ftTOinc  = 12.0; //Length 1.0 in = 0.0254 m  (exact)
const double incTOmm = 25.4; //Length 1.0 in = 25.4 mm  (exact)
const double lbmTOkg = 0.45359237; //Mass 1.0 lbm = 0.45359237 kg (exact)
const double lbmolTOkmol = 0.45359237; //Moles 1.0 lbmol = 0.45359237 kmol (exact)
const double psiTOmpa = 0.006894757; //Pressure 1.0 psia = 0.006894757 MPa
const double mpaTOpsi = 1.0/psiTOmpa;
const double psiTObar = 0.06894757; //Pressure 1.0 psia = 0.006894757 MPa
const double barTOpsi = 1.0/psiTObar;
const double psiTOkpa = 6.894757; //Pressure 1.0 psia = 0.006894757 MPa
const double kpaTOpsi = 1.0/psiTOkpa;

const double atm = 0.101325;  //MPa
const double atmUS = 14.696;  //Psi
const double atmMMHG = 760.0;  //mm Hg
const double atmINCHG = 29.9212;  //inc Hg

double Tcon(char in, double val, char out);

template <typename T>
void bubbleSort(T arr[], int n);

int cubicRoot(double a1, double b1, double c1, double d1, double root[], double& img);
#endif

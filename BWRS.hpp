/**
    This C++ code is the implementation of the analyses presented in the paper
    I.Bell and A. Jäger, "Helmholtz energy translations for common cubic equations of state
    for use in one-fluid and multi-fluid mixture models", J. Res. NIST, 2016

    This code is in the public domain, though if used in academic work, we would appreciate
    a reference back to the paper given above.

 */

//  **********  This code is preliminary, and will be updated.  **********
//  **********  Use only for beta testing.  **********

#ifndef BRW_H
#define BRW_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include "MathUtil.h"

class BWRS
{
protected:
    Kcontribution kcontribution;
    Kconst kconst;
    enum {nSos = 21, Ng = 31}; //for set up as per AGA of ISO
    std::vector<double> Tc, ///< Vector of critical temperatures (in K)
                        rhoc, ///< Vector of critical density
                        w, ///< Vector of acentric factors (unitless)
                        MMi,
                        Ai, //Universal Parameters
                        Bi, //Universal Parameters
                        sosGr; //composition by group of each substance

    std::vector<std::vector<double>> kAi, //Ai data for group contribution
                                     kBi, //Bi data for group contribution
                                     sosGroup; //Contribution of each group for each substance



    double R; ///< The universal gas constant  in J/(mol*K)
    double Mm; //Mm - Molar mass (g/mol)
    size_t N;
    std::vector<double> x; // Molar composition
    //std::string gas[nSos]= {"Methane", "Nitrogen", "Carbon_dioxide", "Ethane", "Propane", "Water", "Hydrogen_sulfide", "Hydrogen", "Carbon_monoxide", "Oxygen",
    //                "iso-Butane", "n-Butane", "iso-Pentane", "n-Pentane", "n-Hexane", "n-Heptane", "n-Octane", "n-Nonane", "n-Decane", "Helium", "Argon"};

    std::vector<std::string> gas = {"74-82-8" ,        "7727-37-9",       "124-38-9",         "74-84-0",      "74-98-6",   "7732-18-5",
                                    "7783-06-4",        "1333-74-0",       "630-08-0",         "7782-44-7",
                                    "75-28-5",          "106-97-8",        "78-78-4",          "109-66-0",     "110-54-3",
                                    "142-82-5",         "111-65-9",        "111-84-2",         "124-18-5",     "7440-59-7", "7440-37-1"};

    double A0;    //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double B0;    //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double C0;    //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double D0;    //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double E0;    //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double a;     //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double b;     //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double c;     //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double d;     //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double ALPHA; //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.
    double GAMMA; //!< A coefficient used when evaluating the BWRS-equation. Independent of pressure and temperature.

    std::vector< std::vector<double> > k;///< The interaction parameters (k_ii = 0)

public:
    BWRS(){

        size_t i,j;
        bool res;
        N = nSos;
        k.resize(N, std::vector<double>(N, 0));
        x.resize(N, 0.0);
        kAi.resize(Ng, std::vector<double>(Ng, 0));
        kBi.resize(Ng, std::vector<double>(Ng, 0));
        sosGroup.resize(N, std::vector<double>(Ng, 0));
        R = 8.3160; // gas constant from Starling [m3 Pa / K mol]

        MMi.resize(N, 0.0);
        Tc.resize(N, 0.0);
        w.resize(N, 0.0);
        rhoc.resize(N, 0.0);

        for(int i= 0; i < N; i++){
            MMi[i]      = std::any_cast<double>(db.get(gas[i], DB::Eq::GEN, DB::DataID::Mm,       res));
            Tc[i]       = std::any_cast<double>(db.get(gas[i], DB::Eq::GEN, DB::DataID::Tc,       res));
            w[i]        = std::any_cast<double>(db.get(gas[i], DB::Eq::GEN, DB::DataID::Acentric, res));
            rhoc[i]     = std::any_cast<double>(db.get(gas[i], DB::Eq::GEN, DB::DataID::Rhoc,       res));
        }

        Ai  = {0.4436900, 1.2843800, 0.3563060, 0.5449790, 0.5286290, 0.4840110, 0.0705233, 0.5040870 , 0.0307452 , 0.0732828, 0.0064500};
        Bi  = {0.115449,  -0.920731, 1.708710, -0.270896,  0.349621,  0.754130, -0.044448,  1.322450,   0.179433,   0.463492, -0.022143};
    };

    double MolarMass(std::map<std::string, double>& mix, bool normalize = true){
        for(std::size_t i = 0; i < N; ++i) x[i]=0.0;
        for(std::size_t i = 0; i < N; ++i) x[i] = mix[gas[i]];

        double xiTot = 0.0;
        for(std::size_t i = 0; i < N; ++i) xiTot += x[i];
        if (normalize && abs(xiTot)>=epsD) for(std::size_t i = 0; i < N; ++i) x[i] /= xiTot;
        Mm = 0;
        for(std::size_t i = 0; i < N; i++) Mm  += x[i]*MMi[i]; // (g*mol-1)
        return Mm;
    }

    /// Get a constant reference to a constant vector of Mathias-Copeman constants
    double GetMW(){ return Mm; }

void Properties(const double T, const double rho, double &P, double &Z, double &dPdD, double &d2PdD2, double &d2PdTD, double &dPdT, double &U, double &H, double &S, double &Cv, double &Cp, double &W, double &G, double &JT, double &Kappa, double &A)
{
    // Sub PropertiesGERG(T, D, x, P, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa, A)

    // Inputs:
    //      T - Temperature (K)
    //      D - Density (g/m^3)

    // Outputs:
    //      P - Pressure (MPa)
    //      Z - Compressibility factor
    //   dPdD - First derivative of pressure with respect to density at constant temperature [kPa/(mol/l)]
    // d2PdD2 - Second derivative of pressure with respect to density at constant temperature [kPa/(mol/l)^2]
    // d2PdTD - Second derivative of pressure with respect to temperature and density [kPa/(mol/l)/K]
    //   dPdT - First derivative of pressure with respect to temperature at constant density (kPa/K)
    //      U - Internal energy (J/mol)
    //      H - Enthalpy (J/mol)
    //      S - Entropy [J/(mol-K)]
    //     Cv - Isochoric heat capacity [J/(mol-K)]
    //     Cp - Isobaric heat capacity [J/(mol-K)]
    //      W - Speed of sound (m/s)
    //      G - Gibbs energy (J/mol)
    //     JT - Joule-Thomson coefficient (K/kPa)
    //  Kappa - Isentropic Exponent
    //      A - Helmholtz energy (J/mol)

    calculateCoefficients(T);

    const double T2 = T*T;
    const double T3 = T2*T;
    const double T4 = T3*T;
    const double T5 = T4*T;
    const double T6 = T5*T;
    const double T7 = T6*T;
    const double T8 = T7*T;

    const double rho2 = rho*rho;
    const double rho3 = rho2*rho;
    const double rho5 = rho2*rho3;
    const double rho6 = rho5*rho;

    const double R2 = R*R;
    const double R3 = R2*R;
    const double R5 = R3*R2;
    const double R6 = R5*R;

    const double RT = R*T;
    const double RT2 = RT*RT;
    const double RT3 = RT2*RT;
    const double RT5 = RT3*RT2;
    const double RT6 = RT5*RT;

    P = calcP(T, rho);

    const double P2 = P*P;
    const double P3 = P2*P;
    const double P4 = P3*P;
    const double P5 = P4*P;

    // Z factor
    Z = P/(rho*R*T);

    // dzdT_P (derivative of Z w.r.t. temperature, at constant pressure)
    const double col = Z;
    const double col2 = col*col;
    const double col3 = col2*col;
    const double col4 = col3*col;
    const double col5 = col4*col;
    const double colgasConstT2 = pow((col*R*T), 2);
    const double expo = exp(-GAMMA*rho2);

    const double nom1 =
          (-B0/(R*T2) + 2*A0/(R2*T3)+ 4*C0/(R2*T5) - 5*D0/(R2*T6) + 6*E0/(R2*T7)) *col4*P
        + (- 2*b/(R2*T3) +  3*a/(R3*T4) +  4*d/(R3*T5)) *col3*P2
        + ALPHA*(-6*a/(R6*T7) - 7*d/(R6*T8)) *P5
        + col3*c*P2/R3 * (-5/T6 - 7*GAMMA*P2/(col2*R2*T8)) * expo
        + col3*c*P2/R3 * ( 1/T5 +   GAMMA*P2/(col2*R2*T7)) * expo*2*GAMMA*P2/(colgasConstT2*T);

    const double den1 =
            6*col5
            - 5*col4
            - ( B0/(R*T) - A0/(R2*T2) - C0/(R2*T4) + D0/(R2*T5) - E0/(R2*T6) )*4*col3*P
            - ( b/(R2*T2) - a/(R3*T3) - d/(R3*T4) )*3*col2*P2
            - 3*col2*c*P2/(R3*T5)*expo
            - GAMMA*P4*c/(R5*T7)*expo
            - col3*c*P2/R3*(1/T5 + GAMMA*P2/(col2*R2*T7))*expo*2*GAMMA*P2/(colgasConstT2*col);
    double dzdT_P = nom1/den1;

    // dzdp_T (derivative of Z w.r.t. pressure, at constant temperature)
    const double nom2 = (B0*R*T - A0 - C0/T2 + D0/T3 - E0/T4) * col4/RT2 + (b*R*T - a - d/T)*col3*2*P/RT3
            + ALPHA*(a + d/T)*5*P4/RT6 + expo*(c*2*P*col3/(T2*RT3) + 4*P3*c*col*GAMMA/(T2*RT5))
            - GAMMA*2*P/colgasConstT2*expo*(c*P2*col3/(T2*RT3) + c*P4*col*GAMMA/(T2*RT5));
    const double den2 = 6*col5 - 5*col4 - (B0*R*T - A0 - C0/T2 + D0/T3 - E0/T4)*4*col3*P/RT2
            - (b*R*T - a - d/T)*3*col2*P2/RT3 - expo*(3*col2*c*P2/(T2*RT3) + c*P4*GAMMA/(T2*RT5))
            - expo*2*GAMMA*P2/(col3*RT2)*(c*P2*col3/(T2*RT3) + c*P4*col*GAMMA/(T2*RT5));
    double dzdp_T = nom2/den2;

    // dzdT_rho (derivative of Z w.r.t. temperature, at constant density)
    double dzdT_rho = A0*rho/(R*T2)
            + 3*C0*rho/(R*T4)
            - 4*D0*rho/(R*T5)
            + 5*E0*rho/(R*T6)
            + a*rho2/(R*T2)
            + 2*d*rho2/(R*T3)
            - ALPHA*a*rho5/(R*T2)
            - 2*ALPHA*d*rho5/(R*T3)
            - 3*c*rho2/(R*T4)*(1 + GAMMA*rho2)*(exp(-GAMMA*rho2));

    Cp = 0.0;
    Cv = 0.0;
    return;
}

//
// T (K)
// P (Mpa.a)
// D (g/m3)
int Density(const int iFlag, const double T, const double P, double &D, int &ierr, std::string &herr){
int I;
double X3, F, F1, F2, F3;
double TOL = 0.5e-9;

    calculateCoefficients(T);
    double X1 = 0.000001, X2 = 1.0; //Density
    D = 0.0;
    F1 = calcP(T, X1);
    F2 = calcP(T, X2);
    F1 -= P;
    F2 -= P;
    if (F1*F2>=0){
            ierr = 1;
            herr = "Root Not inside Interval";
            return 0;
        }

    //------------------
    // BEGIN ITERATING
    //------------------
    for (I=1; I<=50; I++){
    // Use False Position to get point 3.*/
        X3 = X1-F1*(X2-X1)/(F2-F1);
        F3=calcP(T, X3)-P;
    // Use points 1, 2, and 3 to estimate the root using Chamber's method (quadratic solution).
        D = X1*F2*F3/((F1-F2)*(F1-F3)) + X2*F1*F3/((F2-F1)*(F2-F3)) + X3*F1*F2/((F3-F1)*(F3-F2));
        if ((D-X1)*(D-X2)>=0) D=(X1+X2)/2;
        F=calcP(T, D)-P;
        if (fabs(F)<=TOL){
            ierr = 0;
            herr = "";
            return I;
        }
    // Discard quadratic solution if false position root is closer.
        if (fabs(F3) < fabs(F) && F*F3 > 0){
            if (F3*F1>0) {X1=X3; F1=F3;}
            else {X2=X3; F2=F3;}}
        else{
    // Swap in new value from quadratic solution
            if (F*F3 < 0) {X1=D; F1=F; X2=X3; F2=F3;}
            else if (F3*F1 > 0) {X1=D; F1=F;}
            else {X2=D; F2=F;}}
    }
    ierr = 2;
    herr = "Fail to converge";
    D = 0.0;
    return I;
}

protected:

    inline double expW(size_t i){return exp(-3.8*MMi[i]);}

    double calcP(double T, double rho){
    const double T2 = T*T;
    const double T3 = T2*T;
    const double T4 = T3*T;
    const double T5 = T4*T;
    const double T6 = T5*T;
    const double T7 = T6*T;
    const double T8 = T7*T;

    const double rho2 = rho*rho;
    const double rho3 = rho2*rho;
    const double rho5 = rho2*rho3;
    const double rho6 = rho5*rho;

        double P = rho*R*T + (B0*R*T-A0-C0/T2+D0/T3-E0/T4)*rho2+
                   (b*R*T-a-d/T)*rho3 + ALPHA*(a+d/T)*rho6+
                   c*rho3/T2*(1.0+GAMMA*rho2)*exp(-GAMMA*rho2);
        return P;
    }

    double calcdPdD(double T, double rho){
    const double T2 = T*T;
    const double T3 = T2*T;
    const double T4 = T3*T;
    const double T5 = T4*T;
    const double T6 = T5*T;
    const double T7 = T6*T;
    const double T8 = T7*T;

    const double rho2 = rho*rho;
    const double rho3 = rho2*rho;
    const double rho4 = rho3*rho;
    const double rho5 = rho2*rho2*rho;
    const double rho6 = rho5*rho;

        double dPdD = R*T
                    + (B0*R*T - A0 - C0/T2 + D0/T3 - E0/T4)*2.0*rho
                    + (b*R*T - a - d/T)*3.0*rho2
                    + ALPHA*(a + d/T)*6.0*rho5
                    + c*rho2/T2*(3.0 + 3.0*GAMMA*rho2 - 2.0*GAMMA*GAMMA*rho4)*exp(-GAMMA*rho2);
        return dPdD;
    }

    void calculateCoefficients(double T)
    {
        kcontribution.Kij(T,k);
        //kconst.KijUnit(k);
        A0 = 0;
        C0 = 0;
        D0 = 0;
        E0 = 0;
        for (size_t i = 0; i < N; i++)
        {
            double A0i = (Ai[1] + Bi[1]*w[i]) *R       *Tc[i] / rhoc[i];
            double C0i = (Ai[2] + Bi[2]*w[i]) *R  *pow3(Tc[i])/ rhoc[i];
            double D0i = (Ai[8] + Bi[8]*w[i]) *R  *pow4(Tc[i])/ rhoc[i];
            double E0i = (Ai[10] + Bi[10]*w[i]*expW(i)) *R  *pow5(Tc[i])/ rhoc[i];

            double sqrtA0i = sqrt(fabs(A0i))*copysign (1.0, A0i);
            double sqrtC0i = sqrt(fabs(C0i))*copysign (1.0, C0i);
            double sqrtD0i = sqrt(fabs(D0i))*copysign (1.0, D0i);
            double sqrtE0i = sqrt(fabs(E0i))*copysign (1.0, E0i);

            for (size_t j = 0; j < N; j++)
            {
                double A0j = (Ai[1] + Bi[1]*w[j]) *R       *Tc[j] / rhoc[j];
                double C0j = (Ai[2] + Bi[2]*w[j]) *R  *pow3(Tc[j])/ rhoc[j];
                double D0j = (Ai[8] + Bi[8]*w[j]) *R  *pow4(Tc[j])/ rhoc[j];
                double E0j = (Ai[10] + Bi[10]*w[j]*expW(j))*R *pow5(Tc[j])/rhoc[j];

                double sqrtA0j = sqrt(fabs(A0j))*copysign (1.0, A0j);
                double sqrtC0j = sqrt(fabs(C0j))*copysign (1.0, C0j);
                double sqrtD0j = sqrt(fabs(D0j))*copysign (1.0, D0j);
                double sqrtE0j = sqrt(fabs(E0j))*copysign (1.0, E0j);

                A0 += x[i]*x[j]    *(1.0 - k[i][j])*sqrtA0i*sqrtA0j;
                C0 += x[i]*x[j]*pow3(1.0 - k[i][j])*sqrtC0i*sqrtC0j;
                D0 += x[i]*x[j]*pow4(1.0 - k[i][j])*sqrtD0i*sqrtD0j;
                E0 += x[i]*x[j]*pow5(1.0 - k[i][j])*sqrtE0i*sqrtE0j;
            }
        }

        B0 = 0;
        a = 0;
        b = 0;
        c = 0;
        d = 0;
        ALPHA = 0;
        GAMMA = 0;
        for (size_t i =0 ; i < N; i++)
        {
            double B0i =    (Ai[0] + Bi[0]*w[i]) *1.0            /     rhoc[i];
            double ai =     (Ai[5] + Bi[5]*w[i]) *R       *Tc[i] /pow2(rhoc[i]);
            double bi =     (Ai[4] + Bi[4]*w[i]) *1.0            /pow2(rhoc[i]);
            double ci =     (Ai[7] + Bi[7]*w[i]) *R  *pow3(Tc[i])/pow2(rhoc[i]);
            double di =     (Ai[9] + Bi[9]*w[i]) *R  *pow2(Tc[i])/pow2(rhoc[i]);
            double alphai = (Ai[6] + Bi[6]*w[i]) *1.0            /pow3(rhoc[i]);
            double gammai = (Ai[3] + Bi[3]*w[i]) *1.0            /pow2(rhoc[i]);

            B0 += x[i]*B0i;
            a  += x[i]*cbrt(ai);
            b  += x[i]*cbrt(bi);
            c  += x[i]*cbrt(ci);
            d  += x[i]*cbrt(di);

            ALPHA += x[i]*cbrt(alphai);
            GAMMA += x[i]*sqrt(gammai);
        }
        a = pow3(a);
        b = pow3(b);
        c = pow3(c);
        d = pow3(d);
        ALPHA = pow3(ALPHA);
        GAMMA = pow2(GAMMA);
    }

};

#endif

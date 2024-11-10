/**
    This C++ code is the implementation of the analyses presented in the paper
    I.Bell and A. Jäger, "Helmholtz energy translations for common cubic equations of state
    for use in one-fluid and multi-fluid mixture models", J. Res. NIST, 2016

    This code is in the public domain, though if used in academic work, we would appreciate
    a reference back to the paper given above.

 */

//  **********  This code is preliminary, and will be updated.  **********
//  **********  Use only for beta testing.  **********

#ifndef CUBIC_H
#define CUBIC_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include "MathUtil.h"
#include "Kij.hpp"
#include "DB.hpp"

extern DB db;

class AbstractCubic
{
protected:
    //Kconst KijConst;
    Kcontribution kcontribution;
    const double rho_r =1.0, T_r=1.0;
    std::vector<double> Tc, ///< Vector of critical temperatures (in K)
                        Pc, ///< Vector of critical pressures (in Pa)
                        acentric, ///< Vector of acentric factors (unitless)
                        MMi;
    double R_u; ///< The universal gas constant  in J/(mol*K)
    double Mm; //Mm - Molar mass (g/mol)
    double Delta_1, ///< The first cubic constant
           Delta_2; ///< The second cubic constant
    double Pcm, Tcm; // CriticalPressure and Temperature
    double a0=0.0, Da0Dt=0.0, D2a0Dtt=0.0; // Ideal GAs Properties
    std::size_t N; ///< Number of components in the mixture
    std::vector< std::vector<double> > k;///< The interaction parameters (k_ii = 0)
    bool simple_aii; ///< True if the Mathias-Copeman equation for a_ii is not being used
    std::vector<double> C1,C2,C3;///< The Mathias-Copeman coefficients for a_ii
    enum {nSos = 21}; //for set up as per AGA of ISO
    std::vector<double> x; // Molar composition
    //std::string gas[nSos]= {"Methane", "Nitrogen", "Carbon_dioxide", "Ethane", "Propane", "Water", "Hydrogen_sulfide", "Hydrogen", "Carbon_monoxide", "Oxygen",
    //              "iso-Butane", "n-Butane", "iso-Pentane", "n-Pentane", "n-Hexane", "n-Heptane", "n-Octane", "n-Nonane", "n-Decane", "Helium", "Argon"};

    std::vector<std::string> gas = {"74-82-8" ,        "7727-37-9",       "124-38-9",         "74-84-0",      "74-98-6",   "7732-18-5",
                                    "7783-06-4",        "1333-74-0",       "630-08-0",         "7782-44-7",
                                    "75-28-5",          "106-97-8",        "78-78-4",          "109-66-0",     "110-54-3",
                                    "142-82-5",         "111-65-9",        "111-84-2",         "124-18-5",     "7440-59-7", "7440-37-1"};

public:
    /**
     \brief The abstract base clase for the concrete implementations of the cubic equations of state

     This abstract base class describes the structure that must be implemented by concrete implementations
     of the cubic equations of state (SRK, PR, etc.).  The virtual functions must be implemented by the
     derived classes, the remaining functions are generic and are not dependent on the equation of state,
     so long as it has the formulation given in this work.

     */
    AbstractCubic(std::vector<double> Tc,
                  std::vector<double> Pc,
                  std::vector<double> acentric,
                  double R_u,
                  double Delta_1,
                  double Delta_2,
                  std::vector<double> C1 = std::vector<double>(),
                  std::vector<double> C2 = std::vector<double>(),
                  std::vector<double> C3 = std::vector<double>()
                 )
        : Tc(Tc), Pc(Pc), acentric(acentric), R_u(R_u), Delta_1(Delta_1), Delta_2(Delta_2), C1(C1), C2(C2), C3(C3)
        {
            N = static_cast<int>(Tc.size());
            k.resize(N, std::vector<double>(N, 0));
            x.resize(N, 0.0);
            /// If no Mathias-Copeman coefficients are passed in (all empty vectors), use the predictive scheme for m_ii
            simple_aii = (C1.empty() && C2.empty() && C3.empty());
        };

    AbstractCubic(double R_u,
                  double Delta_1,
                  double Delta_2
                 )
        :  R_u(R_u), Delta_1(Delta_1), Delta_2(Delta_2)
        {
        size_t i,j;
            bool res;
            N = nSos;
            k.resize(N, std::vector<double>(N, 0));
            x.resize(N, 0.0);
            MMi.resize(N, 0.0);
            Tc.resize(N, 0.0);
            Pc.resize(N, 0.0);
            acentric.resize(N, 0.0);

            simple_aii = true;
            for(int i= 0; i < N; i++){
                MMi[i]      = std::any_cast<double>(db.get(gas[i], DB::Eq::GEN, DB::DataID::Mm,       res));
                Tc[i]       = std::any_cast<double>(db.get(gas[i], DB::Eq::GEN, DB::DataID::Tc,       res));
                Pc[i]       = std::any_cast<double>(db.get(gas[i], DB::Eq::GEN, DB::DataID::Pc,       res));
                acentric[i] = std::any_cast<double>(db.get(gas[i], DB::Eq::GEN, DB::DataID::Acentric, res));
            }
            //KijConst.Kij(k);
        };

    void MolarMassXi(std::vector<double>& _x, bool normalize = false){
        for(std::size_t i = 0; i < _x.size(); i++) x[i] = _x[i];

        double xiTot = 0.0;
        for(std::size_t i = 0; i < N; i++) xiTot += x[i];
        if (normalize && abs(xiTot)>=epsD) for(std::size_t i = 0; i < N; i++) x[i] /= xiTot;
        Pcm = 0.0; Tcm = 0.0;
        for(std::size_t i = 0; i < N; i++) {
                Pcm += x[i]*Pc[i];
                Tcm += x[i]*Tc[i];
        }
        Pcm *= 1000.0; // Mpa tp Kpa
        return;
    }

    double MolarMass(std::map<std::string, double>& mix, bool normalize = true){
        for(std::size_t i = 0; i < N; ++i) x[i]=0.0;
        for(std::size_t i = 0; i < N; ++i) x[i] = mix[gas[i]];

        double xiTot = 0.0;
        for(std::size_t i = 0; i < N; ++i) xiTot += x[i];
        if (normalize && abs(xiTot)>=epsD) for(std::size_t i = 0; i < N; ++i) x[i] /= xiTot;
        Mm = 0; Pcm = 0.0; Tcm = 0.0;
        for(std::size_t i = 0; i < N; i++) {
                Mm += x[i]*MMi[i]; // (g*mol-1)
                Pcm += x[i]*Pc[i];
                Tcm += x[i]*Tc[i];
        }
        Pcm *= 1000.0; // Mpa tp Kpa
        return Mm;
    }

    /// Copy properties for Ideal Gas
    void SetIdeal(double a0_,double Da0Dt_,double D2a0Dtt_){
        a0 = a0_;
        Da0Dt = Da0Dt_;
        D2a0Dtt = D2a0Dtt_;
    }
    /// Get the vector of critical temperatures (in K)
    std::vector<double> &get_Tc(){ return Tc; }
    /// Get the vector of critical pressures (in Pa)
    std::vector<double> &get_pc(){ return Pc; }
    /// Get the vector of acentric factors
    std::vector<double> &get_acentric(){ return acentric; }
    /// Read-only accessor for value of Delta_1
    double get_Delta_1(){ return Delta_1; }
    /// Read-only accessor for value of Delta_2
    double get_Delta_2(){ return Delta_2; }
    /// Read-only accessor for value of R_u (universal gas constant)
    double get_R_u(){ return R_u; }
    /// Get a constant reference to a constant vector of Mathias-Copeman constants
    double GetMW(){ return Mm; }
    const std::vector<double> &get_C_ref(int i){
        switch (i){
        case 1: return C1;
        case 2: return C2;
        case 3: return C3;
        default:
            throw -1;
        }
    }
    double get_Tr(){ return T_r; }
    double get_rhor(){ return rho_r; }

void Properties(const double T, const double rho, double &P, double &Z, double &dPdD, double &d2PdD2, double &d2PdTD, double &dPdT, double &U, double &H, double &S, double &Cv, double &Cp, double &W, double &G, double &JT, double &Kappa, double &A)
{
    // Sub PropertiesGERG(T, D, x, P, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa, A)

    // Inputs:
    //      T - Temperature (K)
    //      D - Density (Kmol/m^3)
    //    x() - Composition (mole fraction)

    // Outputs:
    //      P - Pressure (kPa)
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
    kcontribution.Kij(T,k);
    double ar, DarDd, D2arDdd, DarDt, D2arDtt, D2arDtDd, phi1, phi2, a, DaDt, D2aDtt;
    double tau = T_r/T, delta = rho/1000.0/rho_r; //from kmol/kg to mol/kg
    ar = alphar(tau, delta, 0, 0);
    DarDd = alphar(tau, delta, 0, 1);
    D2arDdd = alphar(tau, delta, 0, 2);
    DarDt = alphar(tau, delta, 1, 0);
    D2arDtt = alphar(tau, delta, 2, 0);
    D2arDtDd = alphar(tau, delta, 1, 1);

    a = a0+ar;
    DaDt = Da0Dt+tau*DarDt;
    D2aDtt = D2a0Dtt+tau*tau*D2arDtt;

    phi1 = 1.0+2.0*delta*DarDd+delta*delta*D2arDdd;
    phi2 = 1.0+delta*DarDd-tau*delta*D2arDtDd;

    P = rho*R_u*T*(1.0+delta*DarDd);         //Pressure (kPa)
    Z = 1.0+delta*DarDd;         //Compressibility factor
    dPdD = R_u*T*(1.0+2.0*delta*DarDd+delta*delta*D2arDdd);;      //First derivative of pressure with respect to density at constant temperature [kPa/(mol/l)]
    d2PdD2 = 0.0;    //Second derivative of pressure with respect to density at constant temperature [kPa/(mol/l)^2]
    d2PdTD = 0.0;    //Second derivative of pressure with respect to temperature and density [kPa/(mol/l)/K]
    dPdT = rho*R_u*(1.0+delta*DarDd-delta*tau*D2arDtt);      //First derivative of pressure with respect to temperature at constant density (kPa/K)
    U = R_u*T*DaDt;         //Internal energy (J/mol)
    H = R_u*T*(1.0+delta*DarDd+DaDt);        //Enthalpy (J/mol)
    S = R_u*(DaDt-a);         //Entropy [J/(mol-K)]
    Cv = -R_u*D2aDtt;         //Isochoric heat capacity [J/(mol-K)]
    Cp = R_u*(-D2aDtt+phi2*phi2/phi1);         //Isobaric heat capacity [J/(mol-K)]
    W = R_u*T/(Mm/1000.0)*(phi1-phi2*phi2/D2aDtt);         //Speed of sound (m/s)  -  Mm is g/mol we need Kg/mol
    W = sqrt(W);
    G = R_u*T*(1.0+delta*DarDd+a);        //Gibbs energy (J/mol)
    JT = ((phi2-phi1)/(phi2*phi2-D2aDtt*phi1))/(R_u*rho);       //Joule-Thomson coefficient (K/kPa)
    Kappa =phi1/(1.0+delta*DarDd)*(1.0-phi2*phi2/(D2aDtt*phi1));    // Isentropic Exponent
    A = R_u*T*a;        //Helmholtz energy (J/mol)
}

protected: //change to PUBLIC for test program
    /// Get the leading constant in the expression for the pure fluid attractive energy term
    /// (must be implemented by derived classes)
    virtual double a0_ii(std::size_t i) = 0;
    /// Get the leading constant in the expression for the pure fluid covolume term
    /// (must be implemented by derived classes)
    virtual double b0_ii(std::size_t i) = 0;
    /// Get the m_ii variable in the alpha term inculuded in the attractive part
    virtual double m_ii(std::size_t i) = 0;

    /// The residual non-dimensionalized Helmholtz energy \f$\alpha^r\f$
    double alphar(double tau, double delta, std::size_t itau, std::size_t idelta);
    /// The first composition derivative of \f$\alpha^r\f$ as well as derivatives with respect to \f$\tau\f$ and \f$\delta\f$
    double d_alphar_dxi(double tau, double delta,std::size_t itau, std::size_t idelta, std::size_t i, bool xN_independent);
    /// The second composition derivative of \f$\alpha^r\f$ as well as derivatives with respect to \f$\tau\f$ and \f$\delta\f$
    double d2_alphar_dxidxj(double tau, double delta, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);
    /// The third composition derivative of \f$\alpha^r\f$ as well as derivatives with respect to \f$\tau\f$ and \f$\delta\f$
    double d3_alphar_dxidxjdxk(double tau, double delta, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);

    /**
     * \brief The n-th derivative of \f$a_m\f$ with respect to \f$\tau\f$
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param x The vector of mole fractions
     * \param itau The number of derivatives of \f$a_m\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_m, itau=1 is d(a_m)/d(tau), etc.)
     */
    virtual double am_term(double tau, std::size_t itau);
    /**
     * \brief The first composition derivative of \f$a_m\f$ as well as derivatives with respect to \f$\tau\f$
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param x The vector of mole fractions
     * \param itau The number of derivatives of \f$a_m\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_m, itau=1 is d(a_m)/d(tau), etc.)
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d_am_term_dxi(double tau, std::size_t itau, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$a_m\f$ as well as derivatives with respect to \f$\tau\f$
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param x The vector of mole fractions
     * \param itau The number of derivatives of \f$a_m\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_m, itau=1 is d(a_m)/d(tau), etc.)
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d2_am_term_dxidxj(double tau, std::size_t itau, std::size_t i, std::size_t j, bool xN_independent);
    /**
     * \brief The third composition derivative of \f$a_m\f$ as well as derivatives with respect to \f$\tau\f$
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param x The vector of mole fractions
     * \param itau The number of derivatives of \f$a_m\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_m, itau=1 is d(a_m)/d(tau), etc.)
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d3_am_term_dxidxjdxk(double tau, std::size_t itau, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);

    /**
     * \brief The term \f$b_{\rm m}\f$ (mixture co-volume)
     * \param x The vector of mole fractions
     */
    virtual double bm_term();
    /** \brief The first composition derivative of \f$b_m\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d_bm_term_dxi(std::size_t i, bool xN_independent);
    /** \brief The second composition derivative of \f$b_m\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d2_bm_term_dxidxj(std::size_t i, std::size_t j, bool xN_independent);
    /** \brief The third composition derivative of \f$b_m\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    virtual double d3_bm_term_dxidxjdxk(std::size_t i, std::size_t j, std::size_t k, bool xN_independent);

    /**
     * \brief The n-th \f$\tau\f$ derivative of \f$a_{ij}(\tau)\f$
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param i The first index
     * \param j The second index
     * \param itau The number of derivatives of \f$a_{ij}\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_{ij}, itau=1 is d(a_ij)/d(tau), etc.)
     */
    double aij_term(double tau, std::size_t i, std::size_t j, std::size_t itau);
    /** The n-th tau derivative of \f$u(\tau)\f$, the argument of sqrt in the cross aij term
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param i The first index
     * \param j The first index
     * \param itau The number of derivatives of \f$a_{ij}\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_{ij}, itau=1 is d(a_ij)/d(tau), etc.)
     */
    double u_term(double tau, std::size_t i, std::size_t j, std::size_t itau);
    /** Take the n-th tau derivative of the \f$a_{ii}(\tau)\f$ pure fluid contribution
     * \param tau The reciprocal reduced temperature \f$\tau=T_r/T\f$
     * \param i The index of the component
     * \param itau The number of derivatives of \f$u\f$ to take with respect to \f$\tau\f$ (itau=0 is just a_{ij}, itau=1 is d(a_ij)/d(tau), etc.)
     */
    double aii_term(double tau, std::size_t i, std::size_t itau);

    /**
     * \brief The term \f$ \psi^{(-)}\f$ and its \f$\tau\f$ and \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     */
    double psi_minus(double delta, std::size_t itau, std::size_t idelta);
    /**
     * \brief The third composition derivative of \f$ \psi^{(-)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
     double d_psi_minus_dxi(double delta, std::size_t itau, std::size_t idelta, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$ \psi^{(-)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
     double d2_psi_minus_dxidxj(double delta, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);
    /**
     * \brief The third composition derivative of \f$ \psi^{(-)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d3_psi_minus_dxidxjdxk(double delta, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);

    /**
     * \brief The term \f$ \Pi_{12}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     *
     * \f[ \Pi_{12} = (1+\Delta_1\bm\rhor \delta)(1+\Delta_2\bm\rhor \delta) \f]
     */
    double PI_12(double delta, std::size_t idelta);
    /**
     * \brief The first composition derivative of \f$ \Pi_{12}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_PI_12_dxi(double delta, std::size_t idelta, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$ \Pi_{12}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_PI_12_dxidxj(double delta, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);
    /**
     * \brief The third composition derivative of \f$ \Pi_{12}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d3_PI_12_dxidxjdxk(double delta, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);

    /**
     * \brief The term \f$ \tau\cdot a_m(\tau)\f$ and its \f$ \tau \f$ derivatives
     * \param tau The reciprocal reduced temperature \f$\tau = \frac{T_c}{T}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     */
    double tau_times_a(double tau, std::size_t itau);
    /**
     * \brief The first composition derivative of \f$ \tau\cdot a_m(\tau)\f$ and its \f$ \tau \f$ derivatives
     * \param tau The reciprocal reduced temperature \f$\tau = \frac{T_c}{T}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_tau_times_a_dxi(double tau, std::size_t itau, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$ \tau\cdot a_m(\tau)\f$ and its \f$ \tau \f$ derivatives
     * \param tau The reciprocal reduced temperature \f$\tau = \frac{T_c}{T}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_tau_times_a_dxidxj(double tau, std::size_t itau, std::size_t i, std::size_t j, bool xN_independent);
    /**
     * \brief The third composition derivative of \f$ \tau\cdot a_m(\tau)\f$ and its \f$ \tau \f$ derivatives
     * \param tau The reciprocal reduced temperature \f$\tau = \frac{T_c}{T}\f$
     * \param x The vector of mole fractions
     * \param itau How many derivatives to take with respect to \f$\tau\f$
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d3_tau_times_a_dxidxjdxk(double tau, std::size_t itau, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);

    /**
     * \brief The term \f$ \psi^{(+)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     *
     * \f[  \psi^{(+)} = \dfrac{\ln\left(\dfrac{\Delta_1\bm\rhor \delta+1}{\Delta_2\bm\rhor \delta+1}\right)}{\bm(\Delta_1-\Delta_2)}  \f]
     */
    double psi_plus(double delta, std::size_t idelta);
    /**
     * \brief The first composition derivative of \f$ \psi^{(+)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_psi_plus_dxi(double delta, std::size_t idelta, std::size_t i, bool xN_independent);
    /**
     * \brief The second composition derivative of \f$ \psi^{(+)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_psi_plus_dxidxj(double delta, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent);
    /**
     * \brief The third composition derivative of \f$ \psi^{(+)}\f$ and its \f$ \delta \f$ derivatives
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param idelta How many derivatives to take with respect to \f$\delta\f$
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d3_psi_plus_dxidxjdxk(double delta, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent);

    /** \brief The term \f$c\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     *
     * \f$c\f$ is given by
     * \f[
     * c = \frac{1}{b_m}
     * \f]
     * \param x The vector of mole fractions
     */
    double c_term(){
        return 1/bm_term();
    };
    /**
     * \brief The first composition derivative of the term \f$c\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_c_term_dxi(std::size_t i, bool xN_independent){
        return -d_bm_term_dxi(i,xN_independent)/pow(bm_term(), 2);
    };
    /**
     * \brief The second composition derivative of the term \f$c\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_c_term_dxidxj(std::size_t i, std::size_t j, bool xN_independent){
        double b = bm_term();
        return (2*d_bm_term_dxi(i, xN_independent)*d_bm_term_dxi(j, xN_independent) - b*d2_bm_term_dxidxj(i,j,xN_independent))/pow(b, 3);
    };
    /**
     * \brief The third composition derivative of the term \f$c\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d3_c_term_dxidxjdxk(std::size_t i, std::size_t j, std::size_t k, bool xN_independent){
        double b = bm_term();
        return 1/pow(b,4)*(2*b*(d_bm_term_dxi(i, xN_independent)*d2_bm_term_dxidxj(j, k, xN_independent)
                                +d_bm_term_dxi(j, xN_independent)*d2_bm_term_dxidxj(i, k, xN_independent)
                                +d_bm_term_dxi(k, xN_independent)*d2_bm_term_dxidxj(i, j, xN_independent)
                                )
                           - pow(b,2)*d3_bm_term_dxidxjdxk(i,j,k,xN_independent)
                           -6*d_bm_term_dxi(i, xN_independent)*d_bm_term_dxi(j, xN_independent)*d_bm_term_dxi(k, xN_independent)
                           );
    };

    /**
     * \brief The term \f$A\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     *
     * \f[
     * A = \log\left(\frac{\Delta_1\delta\rho_r b_m+1}{\Delta_2\delta\rho_r b+1}\right)
     * \f]
     *
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     */
    double A_term(double delta){
        double b = bm_term();
        return log((Delta_1*delta*rho_r*b+1)/(Delta_2*delta*rho_r*b+1));
    };
    /**
     * \brief The first composition derivative of the term \f$A\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d_A_term_dxi(double delta, std::size_t i, bool xN_independent){
        std::size_t idelta = 0;
        return delta*rho_r*d_bm_term_dxi(i,xN_independent)*(Delta_1-Delta_2)/PI_12(delta, idelta);
    };
    /**
     * \brief The second composition derivative of the term \f$A\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d2_A_term_dxidxj(double delta, std::size_t i, std::size_t j, bool xN_independent){
        std::size_t idelta = 0;
        double PI12 = PI_12(delta, idelta);
        return delta*rho_r*(Delta_1-Delta_2)/pow(PI12, 2)*(PI12*d2_bm_term_dxidxj(i,j,xN_independent)
                                                           - d_PI_12_dxi(delta, 0, j,xN_independent)* d_bm_term_dxi(i,xN_independent));
    };
    /**
     * \brief The third composition derivative of the term \f$A\f$ used in the pure composition partial derivatives of \f$\psi^{(+)}\f$
     * \param delta The reduced density \f$\delta = \frac{\rho}{\rho_c}\f$
     * \param x The vector of mole fractions
     * \param i The first index
     * \param j The second index
     * \param k The third index
     * \param xN_independent True if \f$x_N\f$ is an independent variable, false otherwise (dependent on other \f$N-1\f$ mole fractions)
     */
    double d3_A_term_dxidxjdxk(double delta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent){
        std::size_t idelta = 0;
        double PI12 = PI_12(delta, idelta);
        // The leading factor
        double lead = delta*rho_r*(Delta_1-Delta_2)/pow(PI12, 3);
        return lead*(-PI12*(d_PI_12_dxi(delta, idelta, j, xN_independent)*d2_bm_term_dxidxj(i,k,xN_independent)
                            +d_PI_12_dxi(delta, idelta, k, xN_independent)*d2_bm_term_dxidxj(i,j,xN_independent)
                            +d_bm_term_dxi(i,xN_independent)*d2_PI_12_dxidxj(delta, idelta, j, k, xN_independent))
                    +pow(PI12, 2)*d3_bm_term_dxidxjdxk(i, j, k, xN_independent)
                    +2*d_PI_12_dxi(delta, idelta, j, xN_independent)*d_PI_12_dxi(delta, idelta, k, xN_independent)*d_bm_term_dxi(i, xN_independent)
                    );
    };
};

//==================================================================================================
class PengRobinson : public AbstractCubic
{
public:
    PengRobinson(std::vector<double> Tc,
                 std::vector<double> pc,
                 std::vector<double> acentric,
                 double R_u,
                 std::vector<double> C1 = std::vector<double>(),
                 std::vector<double> C2 = std::vector<double>(),
                 std::vector<double> C3 = std::vector<double>()
                 )
        : AbstractCubic(Tc, pc, acentric, R_u, 1+sqrt(2.0), 1-sqrt(2.0),C1,C2,C3) {};

    PengRobinson(double Tc,
        double pc,
        double acentric,
        double R_u)
        : AbstractCubic(std::vector<double>(1,Tc), std::vector<double>(1,pc), std::vector<double>(1,acentric), R_u, 1+sqrt(2.0), 1-sqrt(2.0)) {};

    PengRobinson()
        : AbstractCubic(8.3144598, 1+sqrt(2.0), 1-sqrt(2.0)) {};


    virtual double a0_ii(std::size_t i);
    virtual double b0_ii(std::size_t i);
    virtual double m_ii(std::size_t i);

    void Density(const int iFlag, const double T, const double P, double &D, int &ierr, std::string &herr){
    // Inputs:
    //  iFlag - Set to 0 for strict pressure solver in the gas phase without checks (fastest mode, but output state may not be stable single phase)
    //          Set to 1 to make checks for possible 2-phase states (result may still not be stable single phase, but many unstable states will be identified)
    //          Set to 2 to search for liquid phase (and make the same checks when iFlag=1)
    //      T - Temperature (K)
    //      P - Pressure (kPa)
    //    x() - Composition (mole fraction)
    // (An initial guess for the density can be sent in D as the negative of the guess for roots that are in the liquid phase instead of using iFlag=2)

    // Outputs:
    //      D - Density (mol/l)
    //          For the liquid phase, an initial value can be sent to the routine to avoid
    //          a solution in the metastable or gas phases.
    //          The initial value should be sent as a negative number.
    //   ierr - Error number (0 indicates no error)
    //   herr - Error message if ierr is not equal to zero

    double A,B,Z;
    double root[3], img;
    int nRoot;
    double tau = T_r/T;

        A = am_term(tau,0)*P/1000.0/pow(R_u*T,2);
        B = bm_term()*P/1000.1/R_u/T;
        Z = 0.0;
        nRoot = cubicRoot(1.0, -1.0+B, A-2.0*B-3.0*B*B, -A*B+B*B+B*B*B, root, img);
        if (nRoot==1) Z = root[0];
        if (nRoot==2) Z = root[1];
        if (nRoot==3) if(iFlag==0) Z = root[0];
        if (nRoot==3) if(iFlag==2) Z = root[2];
        D = P/(R_u*Z*T); // mol/m3
        ierr = 0;
        herr = "";
    }
};

//==================================================================================================
class SRK : public AbstractCubic
{
public:
    SRK(std::vector<double> Tc,
        std::vector<double> pc,
        std::vector<double> acentric,
        double R_u,
        std::vector<double> C1 = std::vector<double>(),
        std::vector<double> C2 = std::vector<double>(),
        std::vector<double> C3 = std::vector<double>())
        : AbstractCubic(Tc, pc, acentric, R_u, 1, 0, C1, C2, C3) {};

    SRK(double Tc,
        double pc,
        double acentric,
        double R_u)
        : AbstractCubic(std::vector<double>(1,Tc), std::vector<double>(1,pc), std::vector<double>(1,acentric), R_u, 1, 0) {};

    SRK()
        : AbstractCubic(8.3144598, 1, 0) {};

    virtual double a0_ii(std::size_t i);
    virtual double b0_ii(std::size_t i);
    virtual double m_ii(std::size_t i);

    void Density(const int iFlag, const double T, const double P, double &D, int &ierr, std::string &herr){
    // Inputs:
    //  iFlag - Set to 0 for strict pressure solver in the gas phase without checks (fastest mode, but output state may not be stable single phase)
    //          Set to 1 to make checks for possible 2-phase states (result may still not be stable single phase, but many unstable states will be identified)
    //          Set to 2 to search for liquid phase (and make the same checks when iFlag=1)
    //      T - Temperature (K)
    //      P - Pressure (kPa)
    //    x() - Composition (mole fraction)
    // (An initial guess for the density can be sent in D as the negative of the guess for roots that are in the liquid phase instead of using iFlag=2)

    // Outputs:
    //      D - Density (mol/l)
    //          For the liquid phase, an initial value can be sent to the routine to avoid
    //          a solution in the metastable or gas phases.
    //          The initial value should be sent as a negative number.
    //   ierr - Error number (0 indicates no error)
    //   herr - Error message if ierr is not equal to zero

    double A,B,Z;
    double root[3], img;
    int nRoot;
    double tau = T_r/T;

        A = am_term(tau,0)*P/1000.0/pow(R_u*T,2);
        B = bm_term()*P/1000.1/R_u/T;

        nRoot = cubicRoot(1.0, -1.0, A-B-B*B, -A*B, root, img);
        if (nRoot==1) Z = root[0];
        if (nRoot==2) Z = root[1];
        if (nRoot==3) if(iFlag==0) Z = root[0];
        if (nRoot==3) if(iFlag==2) Z = root[2];
        D = P/(R_u*Z*T);
        if (nRoot==2){
            ierr=2;
            herr = "WARNING: 2 roots find";
        } else if (nRoot==3){
            ierr=3;
            herr = "WARNING: 3 roots find";

        }else {
            ierr = 0;
            herr = "";
        }
    }
};

#endif

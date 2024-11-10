#ifndef __VALVE_ISA__
#define __VALVE_ISA__

#include <vector>
#include "MathUtil.h"

/* VALVE SIZING ACCORDING TO NSI/ISA-75.01.01-2012 (60534-2-1 MOD)
P1  kPa or bar Inlet absolute static pressure measured at point A
P2  kPa or bar Outlet absolute static pressure measured at point B
DP  kPa or bar P1 - P2 Differential pressure between upstream and downstream pressure taps
Pv  kPa or bar Absolute vapor pressure of the liquid at inlet temperature
Pc  kPa or bar Absolute thermodynamic critical pressure
G   Gas s.g. at S.C.
Cf  Critical Flow Factor
Gf  Gas s.g. at flowing temperature
C   Flow coefficient (Kv, Cv)
Fl  Dimensionless Liquid pressure recovery factor of a control valve without attached fittings
Ff  Dimensionless Liquid critical pressure ratio factor
Qs   m3/h Volumetric flow rate @ Standard Condition
rho1 kg/m3 Density of fluid at P1 and T1
rho0 kg/m3 Density of water at 15°C
T1 K Inlet absolute temperature
D1  mm Inlet Pipe Diameter
D2  mm Outlet Pipe Diameter
D0  mm Orifice Diameter
d   mm Valve Size
dvis m2/s Dinamic Viscosity
*/

class valve{
public:

    valve(){
        init();
    }

    void dataPrintout(std::ofstream& logs){
    std::string outFor;

        verbose(logs, "================================================================================", true);
        verbose(logs, std::format("VALVE SIZING  TAG: {}",tag), true);
        verbose(logs, "================================================================================", true);

        std::cout <<"\033[93m";verbose(logs,"FLUID DATA", true);std::cout <<"\033[0m";
        verbose(logs, std::format("Pv......................... {:<15.4f}  [kPa]     {}\n", Pv, "Absolute vapor pressure of the liquid at inlet temperature"));
        verbose(logs, std::format("rho1....................... {:<15.4f}  [kg/m3]   {}\n", rho1, "Density of fluid at P1 and T1"));
        verbose(logs, std::format("Pc......................... {:<15.4f}  [kPa]     {}\n", Pc, "Absolute thermodynamic critical pressure [kPa]"));
        verbose(logs, std::format("Kinematic vis.............. {:<15.12f}  [(cS]     {}\n", Dvis, "Kinematic viscosity  m2/s (cS) < 1 centistoke = 10–6 m2/s"));
        //verbose(logs, std::format("Tc ........................ {:<15.4f}  [K]       {}\n", Tc, "Absolute thermodynamic critical temperature"));
        verbose(logs, std::format("gamma ..................... {:<15.4f}  [-]       {}\n", gamma, "Specific heat ratio"));
        verbose(logs, std::format("Mol. mass ................. {:<15.4f}  [kg/kmol] {}\n", M, "Molecular mass of flowing fluid"));
        verbose(logs, std::format("Comp. factor............... {:<15.4f}  [-]       {}\n", Z1, "Compressibility factor at inlet conditions"));
        verbose(logs, std::format("Standard comp. factor ..... {:<15.4f}  [-]       {}\n", Zs, "Standard compressibility at 101.325 kPa and 273 K"));
        verbose(logs, std::format("Fgamma .................... {:<15.4f}  [-]       {}\n", Fgamma, "factor Fgamma is used adjust xT"));
        verbose(logs, std::format("Density of fluid std ...... {:<15.4f}  [-]       {}\n", Rhos, "Density at standard conditions"));
        verbose(logs, std::format("Fluid Phase ............... {:<c}                {}\n", phase, "Stream Phase L=Liquid / V=Vapour"));

        std::cout <<"\033[93m";verbose(logs,"PROCESS DATA", true);std::cout <<"\033[0m";
        verbose(logs, std::format("Inlet absolute P .......... {:<15.4f}  [kPa]    {}\n", P1, "Inlet absolute static pressure measured at point A (see Figure 1)"));
        verbose(logs, std::format("Inlet absolute t .......... {:<15.4f}  [K]      {}\n", T1, ""));
        verbose(logs, std::format("Outlet absolute static p .. {:<15.4f}  [kPa]    {}\n", P2, "Outlet absolute static pressure measured at point B (see Figure 1)"));
        verbose(logs, std::format("Ratio of actual pressure .. {:<15.4f}  [-]      {}\n", x,  "(P1-P2)/P1 Ratio of actual delta P to inlet absolute P"));
        verbose(logs, std::format("pressure differential ..... {:<15.4f}  [kPa]    {}\n", DP, "P1-P2 LIQUID, pressure differential"));

        std::cout <<"\033[93m";verbose(logs,"VALVE DATA", true);std::cout <<"\033[0m";
        verbose(logs, std::format("Valve ID .................. {:<3d} \n", valID));
        verbose(logs, std::format("Tag ....................... {} \n", tag));
        verbose(logs, std::format("Type ...................... {} \n", type));
        verbose(logs, std::format("Trim ...................... {} \n", trim));
        verbose(logs, std::format("Flow Direction  ........... {} \n", flowdir));
        verbose(logs, std::format("Nominal valve size ........ {:<15.4f}  [mm]    {}\n", d, "Nominal valve size"));
        verbose(logs, std::format("Internal pipe D upstream .. {:<15.4f}  [mm]    {}\n", D1, "Internal diameter of upstream piping"));
        verbose(logs, std::format("Internal pipe D downstr ... {:<15.4f}  [mm]    {}\n", D2, "Internal diameter of downstream piping"));
        verbose(logs, std::format("Valve style modifier ...... {:<15.4f}  [-]     {}\n", Fd, "Valve style modifier (see Annex A)"));
        verbose(logs, std::format("Pressure diff. ratio factor {:<15.4f}  [-]     {}\n", xt, "Pressure differential ratio factor without attached fittings at choked flow"));
        verbose(logs, std::format("Cmax....................... {:<15.4f}  [m3/h]  {}\n", Cmax, "Multiply CV curve otehrwise put 1.0"));
        verbose(logs,"\n",false);
        verbose(logs,"Stroke  = ",false); for (int i=0; i<nData; i++) verbose(logs, std::format("{:>9.3f} ", strokeTab[i]),false); verbose(logs,"\n",false);
        verbose(logs,"C       = ",false); for (int i=0; i<nData; i++) verbose(logs, std::format("{:>9.3f} ", CTab[i]),false);      verbose(logs,"\n",false);
        verbose(logs,"Fl      = ",false); for (int i=0; i<nData; i++) verbose(logs, std::format("{:>9.3f} ", FlTab[i]),false);     verbose(logs,"\n",false);
    }

    void loadFluidData(double lPv, double lrho1, double lPc, double lDvis, double lgamma, double lM, double lZ1, double lZs, char lphase){
        Pv    = lPv;        // Absolute vapor pressure of the liquid at inlet temperature kPa or bar (psia)
        rho1  = lrho1;      // Density of fluid at P1 and T1  kg/m3 (lbm/ft3)
        Pc    = lPc;        // Absolute thermodynamic critical pressure  kPa or bar (psia)
        Dvis  = lDvis;      // Kinematic viscosity  m2/s (cS) < 1 centistoke = 10–6 m2/s >
        gamma = lgamma;     // Specific heat ratio Dimensionless
        M     = lM;         // Molecular mass of flowing fluid kg/kmol (lbm/lbm-mol)
        Z1    = lZ1;        // Compressibility factor at inlet conditions Dimensionless
        Zs    = lZs;        // Standard compressibility at 101.325 kPa and 273 K Dimensionless
        phase = lphase;     // L or V Size valve for Liquid or Vapour stream
        if (phase == 'L' || phase == 'l') phase = 'L';
        if (phase == 'V' || phase == 'v') phase = 'V';

        // The factor xT is based on air near atmospheric pressure as the flowing fluid with a specific he at
        // ratio of 1.40. If the specific heat ratio for the flowing fluid is not 1.40, the factor Fgamma is used to
        // adjust xT.
        Fgamma = gamma / 1.4;
        // Liquid critical pressure ratio factor, Ff
        // This factor is the ratio of the apparent vena contracta
        // pressure at choked flow conditions to the vapor pressure of the liquid at inlet temperature
        Ff = 0.96-0.28*sqrt(Pv/Pc); // Liquid, Ff is the liquid critical pressure ratio factor.For single component fluids approximated from equation

        Rhos = Ps * M * Zs / (R * Ts);  //Density at standard conditions
    }

    void loadProcData(double lP1, double lT1, double lP2){
        //Qs = lQs;           // Qs Standard volumetric flow rate (see definition 3.2) m3/h (scfh)
                            // Q volumetric flow rate m3/h
        P1 = lP1;           // Inlet absolute static pressure measured at point A (see Figure 1)  kPa or bar (psia)
        T1 = lT1;           // Inlet absolute temperature °K
        P2 = lP2;           // Outlet absolute static pressure measured at point B (see Figure 1) kPa or bar (psia)

        // Caculated Value
        x = (P1-P2) / P1;           // Ratio of actual pressure differential to inlet absolute pressure (P/P1) Dimensionless
        DP = P1-P2;                 // LIQUID, pressure differential
    }

    void  loadValData(int lvalID, std::string ltag, std::string ltype, std::string ltrim, std::string lflowdir, double ld, double lD1, double lD2, double lFd, double lxt, double lCmax, double lCTab[], double lFlTab[]){
    double Csi1, CsiB1, Csi2, CsiB2, Fl;

        valID   = lvalID,
        tag = ltag;
        type = ltype;
        trim = ltrim;
        flowdir = lflowdir;
        d       = ld;       // Nominal valve size mm (in)
        D1      = lD1;      // Internal diameter of upstream piping mm (in)
        D2      = lD2;      // Internal diameter of downstream piping mm (in)
        Fd      = lFd;      // Valve style modifier (see Annex A) Dimensionless
        xt      = lxt;      // Pressure differential ratio factor of a control valve without attached fittings at choked flow Dimensionless
        Cmax    = lCmax;    // Multiply CV curve in [0, 1] otehrwise put 1.0
                            // FlTab Liquid pressure recovery factor of a control valve without attached fittings Dimensionless
                            // CTab  valve C curve multiplied by Cmax
        for (int i=0; i<nData; i++) {CTab[i] = lCTab[i]*Cmax; FlTab[i] = lFlTab[i];}

        if (valID>0)
            for (int i=0; i<valDatas.size(); i++)
                if(valID==valDatas[i].ID){
                    type= valDatas[i].type;
                    trim= valDatas[i].trim;
                    flowdir= valDatas[i].flowDir;
                    xt = valDatas[i].xt;
                    Fd = valDatas[i].Fd;
                    Fl = valDatas[i].Fl;
                    for (int i=0; i<nData; i++) FlTab[i] = Fl;
                }

        // Caculated Value
        Aup     = pi / 4.0 * pow2(D1 / 1000.0);   // Area UPstream m2
        Adown   = pi / 4.0 * pow2(D2 / 1000.0);   // Area DOWNstream m2
        Avalve  = pi / 4.0 * pow2(d / 1000.0);    // Area Valve m2
        Csi1    = 0.5 * pow2(1.0 - pow2(d / D1)); // Upstream velocity head loss coefficient of fitting Dimensionless
        Csi2    = pow2(1.0 - pow2(d / D2));       // Downstream velocity head loss coefficient of fitting Dimensionless
        CsiB1   = 1.0 - pow4(d / D1);             // Inlet Bernoulli coefficient Dimensionless
        CsiB2   = 1.0 - pow4(d / D2);             // Outlet Bernoulli coefficient Dimensionless
        SumCsi  = Csi1 + CsiB1 + Csi2 - CsiB2;
        SumCsi1 = Csi1 + CsiB1;
    }

    //convert from CV to Kv
    void convCvToKv(){
        for(int i =0; i< CTab.size(); i++) CTab[i] *=CvToKv;
    }

    //Restore Stroke from 0 to 0.1 with steps 0.1
    void resetStroke(){
        for (int i=0; i<nData; i++) strokeTab[i] = ((double)i/(nData-1));
    }

    //Get Tabulated Value
    void getTabValue(int ind, double& stroke, double& C, double& Fl){
        stroke = 0.0;
        C = 0.0;
        Fl = 0.0;
        if(ind<0 || ind >= nData) return ;
        stroke = strokeTab[ind];
        C = CTab[ind];
        Fl = FlTab[ind];
    }

    //Restore Stroke from 0 to 0.1 with steps 0.1
    void quarterToLinear(){
    double deltaX, minX, maxX;
    std::vector<double> strokeTabl(10), CTabl(10), FlTabl(10);

        minX = strokeTab[0];
        maxX = strokeTab[nData-2];
        deltaX = maxX-minX;
        // prepare a local copy with range 0 to 0.9
        for(int i = 0; i<nData-1; i++){
            strokeTabl[i] = strokeTab[i];
            CTabl[i] = CTab[i];
            FlTabl[i] = FlTab[i];
        }

        // Resample stroketab
        strokeTab[0] = 0.0;
        for(int i = 1; i<nData; i++) strokeTab[i] = strokeTab[i-1] + deltaX / (double)(nData-1);
        for(int i = 0; i<nData; i++){
            lagrange_value_1d (strokeTab[i], CTab[i],  strokeTabl, CTabl);
            lagrange_value_1d (strokeTab[i], FlTab[i], strokeTabl,  FlTabl);
        }

        //Stroke from 0 to 0.1 with steps 0.1
        for(int i = 0; i<nData; i++)
            strokeTab[i] = (strokeTab[i] - minX) / deltaX;
    }

    void CfromFlow(double& C, double& Qs, double& Q, double& W, double& stroke, double& Re, bool& inScope,double&  vValve){
    double Fp, xtp, Fl, Flp;

        // VAPOUR PHASE
        if (phase == 'V'){
            if (d == D1 && d == D2){
                //Since the valve is line-sized, Fp = 1 and Xtp = Xt.
                Fp = 1.0;
                xtp = xt;
                C = CalcCWOFittComp(Qs, Fp, xtp);
            }else{
                C = CFitting(Q);
            }
            calcQsComp(C, Qs, Q, W);
        } else {

        // LIQUID PHASE
            Qs = 0.0; // NOT USED
            //if (d == D1 && d == D2){
            //    // Since the valve is line-sized, Fp=1 and Flp=Fl.
            //    Fp = 1.0;
            //    Flp = FlTab[0];
            //    C = CalcCWOFittIncomp(Q, Fp, Flp);
            //}else{
                C = CFitting(Q);
            //}
            calcQIncomp(C, Q, W);
        }

        lagrange_value_1d (C, stroke, CTab, strokeTab);
        lagrange_value_1d (C, Fl, CTab, FlTab);
        Re = calcRev(C, Q, Fl);
        inScope = (C / N18/(d*d) < 0.047);
        vValve = Q / 3600.0 / Avalve;
    }

    void Flow(double C, double& Qs, double& Q, double& W, double& stroke, double& Re, bool& inScope,double&  vValve){
    double Fl;

        if (phase == 'V') calcQsComp(C, Qs, Q, W);
        else              calcQIncomp(C, Q, W);

        lagrange_value_1d (C, stroke, CTab, strokeTab);
        lagrange_value_1d (C, Fl, CTab, FlTab);
        Re = calcRev(C, Q, Fl);
        inScope = (C / N18/(d*d) < 0.047);
        vValve = Q / 3600.0 / Avalve;
    }

private:

    struct valData {
        int ID;
        std::string type;
        std::string trim;
        std::string flowDir;
        double Fl;
        double xt;
        double Fd;
    };

    std::vector<valData> valDatas; // default data
    double Pv, rho1, Pc, Dvis, gamma, M, Z1, Zs, Fgamma, Ff, Rhos;
    char phase;
    double P1, P2, T1, x, DP;
    int valID;
    std::string tag, type, trim, flowdir;
    double d, D1, D2, Fd, xt, Cmax;
    double SumCsi, SumCsi1, Aup, Adown, Avalve, Ps, Ts;
    double CvToKv, R;

    double N1,  N2,  N4,  N5,   N6,    N8,  N9,  N90,  N915;
    double N18, N19, N22, N220, N2215, N27, N32, rho0;
    std::vector<double> strokeTab;
    std::vector<double> CTab;
    std::vector<double> FlTab;
    int    nData;

    //
    // LIQUID
    //
    //The fundamental flow model for Incompressible fluids in the turbulent flow regime is given

    void calcQIncomp(double C, double& Q, double& W){
    double DPchoked, DPsizing, Fp, Flp, Fl;

        Fp = calcFp(C);
        lagrange_value_1d (C, Fl, CTab, FlTab);
        Flp = CalcFlp(C, Fl);
        DPchoked = pow2(Flp/Fp)*(P1-Ff*Pv);
        DPsizing = (DP < DPchoked) ? DP : DPchoked;
        Q  = C* N1 * Fp * sqrt(DPsizing/(rho1/rho0)) ;    // m3/h
        W  = Q*rho1;
        return;
    }

    double CalcCWOFittIncomp(double Q, double Fp, double Flp){
        double DPchoked, DPsizing, C;

        DPchoked = pow2(Flp/Fp)*(P1-Ff*Pv);
        DPsizing = (DP < DPchoked) ? DP : DPchoked;
        C = Q / (N1 * Fp) * sqrt((rho1/rho0)/DPsizing);
        return C;
    }

    //
    // VAPOUR
    //
    //The fundamental flow model for compressible fluids in the turbulent flow regime is given
    void calcQsComp(double C, double& Qs, double& Q, double& W){
    double Fp, Xtp, Xchoked, Xsizing, Y;

        Fp = calcFp(C);
        Xtp = CalcXtp(C, Fp);
        Xchoked = Fgamma * Xtp;
        Xsizing = (x < Xchoked) ? x : Xchoked;
        Y = 1.0 - Xsizing / (3.0 * Xchoked);

        Qs = C*N9*Fp*P1*Y*sqrt(Xsizing/(M*T1*Z1));  // standard m3/h
        Q  = Qs * Ps * T1 * Z1 / (P1 * Ts * Zs);    // m3/h
        W  = C*N6*Fp*Y*sqrt(Xsizing*P1*rho1);       //
        return;
    }

    double CalcCWOFittComp(double Qs, double Fp, double Xtp){
        double Xchoked, Xsizing, Y, C;
        double Ff, Flp;

        Xchoked = Fgamma * Xtp;
        Xsizing = (x < Xchoked) ? x : Xchoked;
        Y = 1.0 - Xsizing / (3.0 * Xchoked);
        C = Qs / (N9 * Fp * P1 * Y) * sqrt((M * T1 * Z1)/Xsizing);
        return C;
    }

    //
    // C WITH FITTINGS
    //
    double CFitting(double flow){
    double Rev, Ff, Fp, DPchoked, flowPredicted, C, W, Qs, Q;
    double error;
    int Counter;
    bool inScope;

    double C_lower, Fp_lower, F_lower;
    double C_upper, Fp_upper, F_upper;
    double C_mid,   Fp_mid,   F_mid;

        //A lower flow interval limit per C.2.2 should be set:
        C_lower = 0.0;
        (phase == 'L')? calcQIncomp(C_lower, flowPredicted, W) : calcQsComp(C_lower, flowPredicted, Q, W);
        F_lower = flow - flowPredicted;

        //An upper flow interval limit per C.2.3 should be set:
        if (SumCsi >= 0.0)
            C_upper = 0.075 * d * d * N18;
        else {
            C_upper = 0.99 * d * d * sqrt(-N2 / SumCsi);
            if (C_upper < 0.075 * d * d * N18 ) C_upper = 0.075 * d * d * N18;
        }
        (phase == 'L')? calcQIncomp(C_upper, flowPredicted, W) : calcQsComp(C_upper, flowPredicted, Q, W);

        F_upper = flow - flowPredicted;

        // Check
        if (F_lower * F_upper > 0.0) return 0.0;

        // Iteration
        Counter = 0;

        do{
            C_mid = (C_upper + C_lower) / 2.0;
            (phase == 'L')? calcQIncomp(C_mid, flowPredicted, W) : calcQsComp(C_mid, flowPredicted, Q, W);
            F_mid = flow - flowPredicted;
            //std::cout << Counter<< " "<<printf("%6.2f",C_lower)<< " "<<printf("%6.2f",C_mid)<< " "<<printf("%6.2f",C_upper)<<"\n";
            if (F_mid * F_upper < 0.0) {
                C_lower = C_mid;
                F_lower = F_mid;
            } else {
                C_upper = C_mid;
                F_upper = F_mid;
            }

        } while (abs(C_upper - C_lower) > 0.1 && Counter++ < 60);

        // FAILURE in finding Kv
        if (Counter >= 60) return 0.0;

        // Final Calculation
        C = (C_upper + C_lower) / 2.0;
        (phase == 'L')? calcQIncomp(C, Q, W) : calcQsComp(C, Qs, Q, W);
        return C;
    }

    // The valve Reynolds Number, ReV, is employed to determine whether the flow is fully turbulent.
    // Tests show that flow is fully turbulent when the valve ReV ≥ 10,000.
    // The flow rate is in actual volumetric flow rate units for both incompressible and
    // compressible flows.
    //  The kinematic viscosity, should be evaluated at flow conditions
    double calcRev(double C, double Q, double Fl){
    double Rev;

        if (Dvis * sqrt(C * Fl) == 0.0) return 0.0;
        Rev = N4 * Fd * Q / (Dvis * sqrt(C * Fl));
        Rev *= pow((Fl * Fl * C * C) / (N2 * pow4(d)) + 1.0, 0.25);
        return Rev;
    }

    // LIQUID
    // Estimated combined liquid pressure recovery factor and piping geometry factor with
    // attached fittings, FLP
    double CalcFlp(double C, double Fl){
    double Flpl;

        Flpl = sqrt(1.0 + Fl * Fl / N2 * SumCsi1 * pow2(C / (d * d)));
        Flpl = Fl / Flpl;
        return Flpl;
    }

    // VAPOUR
    // Piping geometry factor Dimensionless
    double calcFp(double C){
    double Fpl, SqrTemp;

        SqrTemp = 1.0 + SumCsi / N2 * pow2(C / (d * d));
        if (SqrTemp < 0.0) return 0.0;
        Fpl = sqrt(SqrTemp);
        Fpl = 1.0 / Fpl;
        return Fpl;
    }

    //  Estimated pressure differential ratio factor with attached fittings, xTP
    // xT is the pressure differential ratio factor of a control valve installed without reducers or other fittings.
    double CalcXtp(double C, double Fp){
    double Xtpl;

        Xtpl = xt / pow2(Fp);
        Xtpl = Xtpl / (1.0 + xt * SumCsi1 / N5 * pow2(C / (d * d)));
        return Xtpl;
    }


    void init(){
        CvToKv = 1.0 / 1.15622; //Convert flow coefficient in US gpm.psi to flow factor in (m³/hr).bar
        N1 = 0.1;
        N2 = 0.0016;
        N4 = 0.0707;
        N5 = 0.0018;
        N6 = 3.16;
        N8 = 1.1;
        N90 = 24.6;
        N915 = 26.0;
        N18 = 0.865;
        N19 = 2.5;
        N220 = 17.3;
        N2215 = 18.4;
        N27 = 0.775;
        N32 = 140;
        rho0 = 999.84; //999.1026
        // Sample Data
        strokeTab.resize(11);
        CTab.resize(11);
        FlTab.resize(11);
        //          (0) stroke       (1) Cv            (2) Fl
        strokeTab[0]  = 0.0;  CTab[0]  = 0.0;  FlTab[0]  = 0.85;
        strokeTab[1]  = 0.1;  CTab[1]  = 0.1;  FlTab[1]  = 0.85;
        strokeTab[2]  = 0.2;  CTab[2]  = 0.2;  FlTab[2]  = 0.84;
        strokeTab[3]  = 0.3;  CTab[3]  = 0.3;  FlTab[3]  = 0.79;
        strokeTab[4]  = 0.4;  CTab[4]  = 0.4;  FlTab[4]  = 0.75;
        strokeTab[5]  = 0.5;  CTab[5]  = 0.5;  FlTab[5]  = 0.71;
        strokeTab[6]  = 0.6;  CTab[6]  = 0.6;  FlTab[6]  = 0.63;
        strokeTab[7]  = 0.7;  CTab[7]  = 0.7;  FlTab[7]  = 0.58;
        strokeTab[8]  = 0.8;  CTab[8]  = 0.8;  FlTab[8]  = 0.56;
        strokeTab[9]  = 0.9;  CTab[9]  = 0.9;  FlTab[9]  = 0.54;
        strokeTab[10] = 1.0;  CTab[10] = 1.0;  FlTab[10] = 0.54;
        Cmax = 521.0;

        nData = 11;
        for (int i=1; i<nData; i++) CTab[i] *= CvToKv; // Cv to Kv

        // Define Standard Condition @ 0 Degree
        Ps = 101.325; //KPa
        Ts = 273.0;  //K
        N9 = N90;
        N22 = N220;

        //Ps = 101.325; //KPa
        //Ts = 288; //K
        //N9 = N915;
        //N22 = N2215;

        //
        valDatas = { {0,"User Input","User Input","User Input",0,0,0} ,
         {10,"Globe, single port","3 V-port plug","Open or close",0.9,0.7,0.48} ,
         {11,"Globe, single port","4 V-port plug","Open or close",0.9,0.7,0.41} ,
         {12,"Globe, single port","6 V-port plug","Open or close",0.9,0.7,0.3} ,
         {13,"Globe, single port","Contoured plug (linear and equal percentage)","Open",0.9,0.72,0.46} ,
         {14,"Globe, single port","Contoured plug (linear and equal percentage)","Close",0.8,0.55,1} ,
         {15,"Globe, single port","60 equal diameter hole drilled cage","Outward or inward",0.9,0.68,0.13} ,
         {16,"Globe, single port","120 equal diameter hole drilled cage","Outward or inward",0.9,0.68,0.09} ,
         {17,"Globe, single port","Characterized cage, 4-port","Outward",0.9,0.75,0.41} ,
         {18,"Globe, single port","Characterized cage, 4-port","Inward",0.85,0.7,0.41} ,
         {100,"Globe, double port","Ported plug","Inlet between seats",0.9,0.75,0.28} ,
         {101,"Globe, double port","Contoured plug","Either direction",0.85,0.7,0.32} ,
         {200,"Globe, angle","Contoured plug (linear and equal percentage)","Open",0.9,0.72,0.46} ,
         {201,"Globe, angle","Contoured plug (linear and equal percentage)","Close",0.8,0.65,1} ,
         {202,"Globe, angle","Characterized cage, 4-port","Outward c)",0.9,0.65,0.41} ,
         {203,"Globe, angle","Characterized cage, 4-port","Inward c) ",0.85,0.6,0.41} ,
         {204,"Globe, angle","Venturi","Close",0.5,0.2,1} ,
         {300,"Globe, small flow trim","V-notch","Open",0.98,0.84,0.7} ,
         {301,"Globe, small flow trim","Flat seat (short travel)","Close",0.85,0.7,0.3} ,
         {302,"Globe, small flow trim","Tapered needle","Open",0.95,0.84,0} ,
         {400,"Rotary","Eccentric spherical plug","Open",0.85,0.6,0.42} ,
         {401,"Rotary","Eccentric spherical plug","Close",0.68,0.4,0.42} ,
         {402,"Rotary","Eccentric conical plug","Open",0.77,0.54,0.44} ,
         {403,"Rotary","Eccentric conical plug","Close",0.79,0.55,0.44} ,
         {500,"Butterfly (centered shaft)","Swing-through (70°)","Either",0.62,0.35,0.57} ,
         {501,"Butterfly (centered shaft)","Swing-through (60°)","Either",0.7,0.42,0.5} ,
         {502,"Butterfly (centered shaft)","Fluted vane (70°)","Either",0.67,0.38,0.3} ,
         {503,"Butterfly (centered shaft)","Offset seat (70°)","Either",0.67,0.35,0.57} ,
         {600,"Ball","Full bore (70°)","Either",0.74,0.42,0.99} ,
         {601,"Ball","Segmented ball","Either",0.6,0.3,0.98} ,
         {700,"Globe and angle","Multistage Multipath 2","Either",0.97,0.812,0} ,
         {701,"Globe and angle","Multistage Multipath 3","Either",0.99,0.888,0} ,
         {702,"Globe and angle","Multistage Multipath 4","Either",0.99,0.925,0} ,
         {703,"Globe and angle","Multistage Multipath 5","Either",0.99,0.95,0} ,
         {704,"Globe and angle","Multistage Single path 2","Either",0.97,0.896,0} ,
         {705,"Globe and angle","Multistage Single path 3","Either",0.99,0.935,0} ,
         {706,"Globe and angle","Multistage Single path 4","Either",0.99,0.96,0} };
    }
};

#endif // VALVE_ISA

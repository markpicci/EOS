#include "conIO.h"
#include "Valve.hpp"

//
//=================================================================================================
// VALVE
//
void sizeValve(valve& v, std::ofstream& logs){
COORD curPos;
char charIn;

std::string Tag, Type, Trim, FlowDir;
double Pv=0.0, rho1=0.0, Pc=0.0, Tc=0.0, Dvis=0.0, gamma=0.0, M=0.0, Z1=0.0, Zs=0.0;
char phase;
double P1=0.0, T1=0.0, P2=0.0;
double valID=0, d=0.0, D1=0.0, D2=0.0, Fd=0.0, xt=0.0, Cmax=0.0, CTab[11], FlTab[11], FlDef;
double C, Qs, Q, W, stroke, Re, Qval, Qsval, vValve;
bool inScope;
int nData = 11;

    //========================================================================================
    // TEST 1
    Pv =70.1; rho1=965.4; Pc=22120.0; Tc=0.0; Dvis=3.26e-7; gamma=0.0; M= 0.0; Z1=0.0; Zs=0.0; phase='L';
    P1=680.0; T1=363.0; P2=220.0;
    //lFlTab[0] = =0.9; 0.9; 0.9; 0.9; 0.9; 0.9; 0.9; 0.9; 0.9; 0.9; 0.9};
    Tag="V-001", Type="globe", Trim="parabolic plug", FlowDir="flow-to-open";
    valID=-1; d=150.0; D1=150.0; D2=150.0; Fd=0.46; xt=0.0; Cmax= 520; FlDef=0.9;
    Qval=360.0;
    Qsval=0.0;
    // TEST 1 END

    // TEST 2
    Pv =70.1; rho1=965.4; Pc=22120.0; Tc=0.0; Dvis=3.26e-7; gamma=0.0; M= 0.0; Z1=0.0; Zs=0.0; phase='L';
    P1=680.0; T1=363.0; P2=220.0;
    Tag="V-002", Type="ball valve", Trim="segmented ball", FlowDir="flow-to-open";
    valID=-1; d=150.0; D1=150.0; D2=150.0; Fd=0.98; xt=0.0; Cmax= 520; FlDef=0.6;
    Qval=360.0;
    Qsval=0.0;
    // TEST 2 END

    // TEST 3
    Pv =0.0; rho1=8.389; Pc=7.387; Tc=304; Dvis=2.526e-6; gamma=1.3; M= 44.01; Z1=0.991; Zs=0.994; phase='V';
    P1=680.0; T1=433.0; P2=450.0;
    Tag="V-003", Type="rotary", Trim="eccentric spherical plug", FlowDir="flow-to-open";
    valID=-1; d=100.0; D1=100.0; D2=100.0; Fd=0.42; xt=0.6; Cmax= 520; FlDef=0.85;
    Qval=0.0;
    Qsval=3800.0;
    // TEST 3 END

    // TEST 4
    Pv =0.0; rho1=8.389; Pc=7.387; Tc=304; Dvis=2.526e-6; gamma=1.3; M= 44.01; Z1=0.991; Zs=0.994; phase='V';
    P1=680.0; T1=433.0; P2=250.0;
    Tag="V-003", Type="rotary", Trim="eccentric spherical plug", FlowDir="flow-to-open";
    valID=-1; d=100.0; D1=100.0; D2=100.0; Fd=0.42; xt=0.6; Cmax= 520; FlDef=0.85;
    Qval=0.0;
    Qsval=3800.0;
    // TEST 4 END

    // TEST 5
    Pv =4.0; rho1=780.0; Pc=22120.0; Tc=0.0; Dvis=3.26e-7; gamma=0.0; M= 0.0; Z1=0.0; Zs=0.0; phase='L';
    P1=3550.0; T1=0.0; P2=2240.0;
    Tag="V-005", Type="Butterfly", Trim="-", FlowDir="-";
    valID=-1; d=101.6; D1=154.1; D2=202.7; Fd=0.46; xt=0.0; Cmax= 1.0; FlDef=-1.0;
    CTab[0]  = 0.0;   CTab[1] = 17.2;  CTab[2] = 50.2;  CTab[3] = 87.8; CTab[4]  = 146.0; CTab[5] =  206.0; CTab[6] = 285.0; CTab[7] = 365.0; CTab[8] = 465.0; CTab[9] = 521.0; CTab[10] =  521.0;
    FlTab[0] = 0.85; FlTab[1] = 0.85; FlTab[2] = 0.84; FlTab[3] = 0.79; FlTab[4] = 0.75;  FlTab[5] = 0.71;  FlTab[6] = 0.63; FlTab[7] = 0.58; FlTab[8] = 0.56; FlTab[9] = 0.54; FlTab[10] = 0.54;
    Qval=750.0;
    Qsval=0.0;
    // TEST 5 END
    //========================================================================================

ReviseFullData:
    clearScreen();
    std::cout <<"================================================================================";std::cout << "\n";

    std::cout <<"\033[93mFLUID DATA\033[0m\n";
    curPos = getCurPos();
        do{
            setCurPos(curPos);
            std::cout <<"Absolute vapor pressure             [KPa.a]   =  "; Pv = get_double("LhGreenF", Pv); std::cout << "\n";
            std::cout <<"Density of fluid at 1               [kg/m3]   =  "; rho1 = get_double("LhGreenF", rho1); std::cout << "\n";
            std::cout <<"Absolute critical pressure          [KPa.a]   =  "; Pc = get_double("LhGreenF", Pc); std::cout << "\n";
            std::cout <<"Absolute critical temperature       [K]       =  "; Tc = get_double("LhGreenF", Tc); std::cout << "\n";
            std::cout <<"Kinematic viscosity                 [m2/s]    =  "; Dvis = get_double("LhGreenF", Dvis); std::cout << "\n";
            std::cout <<"Specific heat ratio                 [-]       =  "; gamma = get_double("LhGreenF", gamma); std::cout << "\n";
            std::cout <<"Molecular mass                      [kg/kmol] =  "; M = get_double("LhGreenF", M); std::cout << "\n";
            std::cout <<"Compressibility factor 1            [-]       =  "; Z1 = get_double("LhGreenF", Z1); std::cout << "\n";
            std::cout <<"Compressibility factor Standard     [-]       =  "; Zs = get_double("LhGreenF", Zs); std::cout << "\n";
            std::cout <<"Stream Phase (L/V)                  [-]       =  "; phase = choice("Liquid (L) or Vapour (V) ","LV",phase);
            charIn=choice("Confirm (y) or Revise (r) ","yr");
        } while(charIn != 'y');
        v.loadFluidData  (Pv, rho1, Pc, Dvis, gamma, M, Z1, Zs, phase);

        std::cout <<"\033[93mPROCESS DATA\033[0m\n";
        curPos = getCurPos();
        do{
            setCurPos(curPos);
            std::cout <<"Inlet absolute static pressure  [kPa.a]   =  "; P1 = get_double("LhGreenF", P1); std::cout << "\n";
            std::cout <<"Inlet absolute temperature      [K]       =  "; T1 = get_double("LhGreenF", T1); std::cout << "\n";
            std::cout <<"Outlet absolute static pressure [kPa.a]   =  "; P2 = get_double("LhGreenF", P2); std::cout << "\n";
            std::cout <<"FLow Rate Standard Condition    [sm3/h]   =  "; Qsval = get_double("LhGreenF", Qsval); std::cout << "\n";
            std::cout <<"FLow Rate                       [m3/h]    =  "; Qval = get_double("LhGreenF", Qval); std::cout << "\n";
            charIn=choice("Confirm (y) or Revise (r) ","yr");
        } while(charIn != 'y');
        v.loadProcData  (P1, T1, P2);

        clearScreen();
        std::cout <<"\033[93mVALVE DATA\033[0m\n";
        curPos = getCurPos();
        do{
            setCurPos(curPos);
            std::cout <<"Valve ID                                     =  "; valID = get_int("LhGreenF", valID); std::cout << "\n";
            std::cout <<"Tag                                          =  "; Tag = get_string("LhGreenF", Tag); std::cout << "\n";
            std::cout <<"Type                                         =  "; Type = get_string("LhGreenF", Type); std::cout << "\n";
            std::cout <<"Trim                                         =  "; Trim = get_string("LhGreenF", Trim); std::cout << "\n";
            std::cout <<"Flow Direction                               =  "; FlowDir = get_string("LhGreenF", FlowDir); std::cout << "\n";

            std::cout <<"Nominal valve size                    [mm]   =  "; d = get_double("LhGreenF", d); std::cout << "\n";
            std::cout <<"Internal diameter pipe upstream       [mm]   =  "; D1 = get_double("LhGreenF", D1); std::cout << "\n";
            std::cout <<"Internal diameter pipe downstr        [mm]   =  "; D2 = get_double("LhGreenF", D2); std::cout << "\n";
            std::cout <<"Valve style modifier Fd               [-]    =  "; Fd = get_double("LhGreenF", Fd); std::cout << "\n";
            std::cout <<"Pressure differential ratio factor xt [-]    =  "; xt = get_double("LhGreenF", xt); std::cout << "\n";
            std::cout <<"KVmax (Cmax)                          [m3/h] =  "; Cmax = get_double("LhGreenF", Cmax); std::cout << "\n";
            std::cout <<"Fl                                    [-]    =  "; FlDef = get_double("LhGreenF", FlDef); std::cout << "\n";
            //for (int i=0; i<nData; i++) CTab[i] = ((double)i/10.0);
            for (int i=0; i<nData; i++){
                std::cout <<"Opening "<<std::format("{:<4.2f}  -  ", ((double)i/10.0))<<" C (KV)              [-]    =  ";
                CTab[i] = get_double("LhGreenF", CTab[i]);
                //std::cout <<"fl "; FlTab[i] = get_double("LhGreenF", FlTab[i]);
                std::cout << "\n";
            }

            if (FlDef>=0.0) for (int i=0; i<nData; i++) FlTab[i] = FlDef;
            charIn=choice("Confirm (y) or Revise (r) ","yr");
        } while(charIn != 'y');
        v.loadValData(-1, Tag, Type, Trim, FlowDir, d, D1, D2, Fd, xt,  Cmax, CTab, FlTab);

        clearScreen();
        v.resetStroke();
        v.quarterToLinear();
        v.dataPrintout(logs);
        Qs = Qsval; Q = Qval; // Qsval and Qval are input by User, Qs and Q are calculated value
        v.CfromFlow( C, Qs, Q, W, stroke, Re, inScope, vValve);

        std::cout <<"\033[93m";verbose(logs,"\nRESULTS", true);std::cout <<"\033[0m";
        verbose(logs, std::format("Kv.............. {:<15.2f} {}\n", C,      "[m3/h]"));
        verbose(logs, std::format("Qs.............. {:<15.4f} {}\n", Qs,     "[sm3/h]"));
        verbose(logs, std::format("Q............... {:<15.4f} {}\n", Q,      "[m3/h]"));
        verbose(logs, std::format("W............... {:<15.4f} {}\n", W,      "[kg/h]"));
        verbose(logs, std::format("stroke.......... {:<15.4f} {}\n", stroke, "[%]"));
        verbose(logs, std::format("Re.............. {:<15.0f} {}\n", Re,     "[-]"));
        verbose(logs, std::format("In Scope........ {}            {}\n", inScope,"[-]"));
        verbose(logs, std::format("Speed in Valve.. {:<15.4f} {}\n", vValve,     "[m/s]"));

        std::cout <<"\033[93m";verbose(logs,"\nINSTALLED FLOW", true);std::cout <<"\033[0m";
        verbose(logs, std::format("Stroke [%]   Kv [m3/h]       Qs [sm3/h]      Q [m3/h]      Gain [dQ%/ds%]    W [kg/h]        speed [m/s]     re [-]       In Scope\n"));
        double Qold = 0.0, Qmax, dQ=0.0, strokeOld =0.0, strokeDelta;
        double strokeTemp,CTemp, FlTemp;

        v.getTabValue (nData-1, strokeTemp, CTemp, FlTemp);
        v.Flow(CTemp, Qs, Q, W, strokeTemp, Re, inScope, vValve);
        Qmax = (Qs>0.0)? Qs: Q;

        for (int i=0; i<nData; i++){
                v.getTabValue (i, strokeTemp,CTemp, FlTemp);
                v.Flow(CTemp, Qs, Q, W, stroke, Re, inScope, vValve);
                strokeDelta = stroke - strokeOld;
                verbose(logs, std::format("{:<11.4f}  {:<15.4f} {:<15.4f} {:<15.4f} {:<15.4f} {:<15.4f} {:<15.4f} {:<12.0f} {} \n", stroke, CTemp, Qs, Q, dQ/strokeDelta, W, vValve , Re, inScope));
                dQ   = (((Qs>0.0)? Qs: Q)-Qold)/Qmax;
                Qold = (Qs>0.0)? Qs: Q;
                strokeOld = stroke;
        }

        charIn=choice("Exit (y) or Revise (r) ","yr");
        if (charIn == 'r') goto ReviseFullData;
}

//
//=================================================================================================
// MAIN
//
void equipmentMain(std::ofstream& logs){
valve v;
char charIn;
std::vector<std::string> menuList = {"(v)alve", "(e)xit"};

    do{
        charIn = conMenu("MAIN MENU",menuList);
        std::cout <<"\n";

        if(charIn=='v') sizeValve(v, logs);

    } while (charIn != 'e');
}

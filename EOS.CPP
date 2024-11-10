#include "mainh.hpp"
#include "conIO.h"
#include "ISO-20765-5.hpp"
#include "ISO-12213-2.hpp"
#include "AGA8 Gross.hpp"
#include "AGA8 Detail.hpp"
#include "GERG2008.hpp"
#include "generalizedcubic.hpp"
#include "IF97.hpp"
#include "BWRS.hpp"

template<typename... Args>
void logEOSdata(std::ofstream& logs, std::string_view fmt_str, Args&&... args){
    std::string outFor = std::vformat(fmt_str, std::make_format_args(args...));
    std::cout << outFor;
    logs << outFor;
}

void EOS(std::ofstream& logs, const size_t activeMix, std::vector<mixSt>& mix) {
double TT, PP, BMIX;
double D, P, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa, A;
double Tin, Pin, MW;
char charIn;
int ierr, dFlag;
std::string herr;
std::string fname;

ISO122132 Aga8;
ISO207655 Visosity;
AGA8Gross  AGA8Gr;
AGA8Detail AGA8Det;
GERG2008   Gerg;
PengRobinson PR;
BWRS Bwrs;

    clearScreen();
    std::vector<std::string> menuList = {"(1) ISO-12213-2","(2) AGA8 Detail","(3) AGA8Gross","(4) GERG 2008","(5) PENG ROBINSON","(6) BWRS","(9) IAPWS R7-97 [2012]","(P) Pressure", "(r)ange", "(e)xit"};

//===================================== MANUAL
    do{
        charIn = conMenu("EOS MENU",menuList);
        std::cout <<"\n";
        bool MixValid = (activeMix >= 0 &&  activeMix < mix.size());
        if(charIn=='1' && MixValid){
            std::cout <<"================================================================================";std::cout << "\n";
            std::cout <<"\033[93mISO-12213-2\033[0m\n";
            std::cout <<"Pressure [MPa.g] = "; Pin = get_double("LhGreenF"); std::cout << "\n";
            std::cout <<"Temperature [C]  = "; Tin = get_double("LhGreenF"); std::cout << "\n";
            logEOSin(logs, "ISO-12213-2", Pin, Tin, atm);
            TT=Tin+273.15;
            PP=Pin+atm; //PRESSURE is ABSOLUTE
            Aga8.MolarMass(mix[activeMix].mix);
            Aga8.DZOFPT(&PP, &TT, &D, &Z, &BMIX); // Density Kmol/m3
            MW=Aga8.GetMW();
            std::cout <<"\033[93mResults\033[0m\n";
            logEOSdata(logs,"Molar Mass......... {:<15.4f} {}\n", MW  ,"[g/mol]");
            logEOSdata(logs,"Density............ {:<15.3f} {}\n", D*MW,"[kg/m3]");
            logEOSdata(logs,"Z.................. {:<15.6f} {}\n", Z   ,"[-]");
            logEOSdata(logs,"Viscosity.......... {:<15.6f} {}\n\n", Visosity.VIS(TT, D, mix[activeMix].mix),"[mPa*s]");
            choice("Exit (y) ","yY");
        }

        if(charIn=='2' && MixValid){
            std::cout <<"================================================================================";std::cout << "\n";
            std::cout <<"\033[93mAGA 8 DETAIL METHOD\033[0m\n";
            std::cout <<"Pressure [MPa.g] = "; Pin = get_double("LhGreenF"); std::cout << "\n";
            std::cout <<"Temperature [C]  = "; Tin = get_double("LhGreenF"); std::cout << "\n";
            logEOSin(logs, "AGA 8 DETAIL METHOD", Pin, Tin, atm);
            TT=Tin+273.15;
            PP=Pin+atm; //PRESSURE is ABSOLUTE
            AGA8Det.MolarMass(mix[activeMix].mix);
            //printf("%16s xi= %6.4f MWi= %9.5f\n", Aga8.GetName(i).c_str(), ex1[6][i],Aga8.GetMW(i));
            MW = AGA8Det.GetMW();
            D = -1.0;
            AGA8Det.Density(TT, PP*1000.0 ,D, ierr, herr);
            std::cout <<"\033[93mCheck\033[0m\n";
            logEOSdata(logs,"Density Error ID... {:<d}\n", ierr);
            logEOSdata(logs,"Density Error...... {:<s}\n", herr);
            AGA8Det.Properties(TT, D, PP, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa,A);
            std::cout <<"\033[93mResults\033[0m\n";
            logEOSdata(logs,"Molar Mass......... {:<15.4f} {}\n", MW,"[g/mol]");
            logEOSdata(logs,"Density............ {:<15.3f} {}\n", D*MW ,"[kg/m3]");
            logEOSdata(logs,"Pressure........... {:<15.4f} {:10s}       :  {:<15.4f} {}\n", PP/1000.0,"[Mpa.a]", PP*kpaTOpsi, "[psi.a]" );
            logEOSdata(logs,"Z.................. {:<15.6f} {}\n", Z ,"[-]");
            logEOSdata(logs,"U.................. {:<15.2f} {}\n", U/MW, "[KJ/Kg]");
            logEOSdata(logs,"H.................. {:<15.2f} {}\n", H/MW, "[KJ/Kg]");
            logEOSdata(logs,"S.................. {:<15.4f} {}\n", S/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cv................. {:<15.4f} {}\n", Cv/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cp................. {:<15.4f} {}\n", Cp/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Sound Speed........ {:<15.2f} {}\n", W, "[m/S]");
            logEOSdata(logs,"G.................. {:<15.2f} {}\n", G/MW, "[KJ/Kg]");
            logEOSdata(logs,"Joule-Thomson...... {:<15.4f} {}\n", JT*1000.0, "[K/MPa]");
            logEOSdata(logs,"Isentropic Exponent {:<15.7f} {}\n", Kappa, "[-]");
            logEOSdata(logs,"Helmholtz.......... {:<15.4f} {}\n\n", A/MW, "[KJ/kPa]");
            choice("Exit (y) ","yY");
        }

        if(charIn=='P' && MixValid){
            std::cout <<"================================================================================";std::cout << "\n";
            std::cout <<"\033[93mAGA 8 DETAIL METHOD\033[0m\n";
            std::cout <<"Density [Kg/m3] = "; D = get_double("LhGreenF"); std::cout << "\n";
            std::cout <<"Temperature [C]  = "; Tin = get_double("LhGreenF"); std::cout << "\n";
            logEOSin(logs, "AGA 8 DETAIL METHOD", Pin, Tin, atm);
            TT=Tin+273.15;
            AGA8Det.MolarMass(mix[activeMix].mix);
            MW = AGA8Det.GetMW(); //[Kg/Kmol]
            D = D/MW; //[Kg/m3] -> [kmol/m3]
            AGA8Det.Properties(TT, D, PP, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa,A);
            std::cout <<"\033[93mResults\033[0m\n";
            logEOSdata(logs,"Molar Mass......... {:<15.4f} {}\n", MW,"[g/mol]");
            logEOSdata(logs,"Density............ {:<15.3f} {}\n", D*MW ,"[kg/m3]");
            logEOSdata(logs,"Pressure........... {:<15.4f} {:10s}       :  {:<15.4f} {}\n", PP/1000.0,"[Mpa.a]", PP*kpaTOpsi, "[psi.a]" );
            logEOSdata(logs,"Z.................. {:<15.6f} {}\n", Z ,"[-]");
            logEOSdata(logs,"U.................. {:<15.2f} {}\n", U/MW, "[KJ/Kg]");
            logEOSdata(logs,"H.................. {:<15.2f} {}\n", H/MW, "[KJ/Kg]");
            logEOSdata(logs,"S.................. {:<15.4f} {}\n", S/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cv................. {:<15.4f} {}\n", Cv/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cp................. {:<15.4f} {}\n", Cp/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Sound Speed........ {:<15.2f} {}\n", W, "[m/S]");
            logEOSdata(logs,"G.................. {:<15.2f} {}\n", G/MW, "[KJ/Kg]");
            logEOSdata(logs,"Joule-Thomson...... {:<15.4f} {}\n", JT*1000.0, "[K/MPa]");
            logEOSdata(logs,"Isentropic Exponent {:<15.7f} {}\n", Kappa, "[-]");
            logEOSdata(logs,"Helmholtz.......... {:<15.4f} {}\n\n", A/MW, "[KJ/kPa]");
            choice("Exit (y) ","yY");
        }

        if(charIn=='3' && MixValid){
            double Hv, Gr, HN, HCH;
            Hv = 0.0;
            double Th; //     Th - Reference temperature for heating value (K)
            double Td; //     Td - Reference temperature for density (K)
            double Pd; //     Pd - Reference pressure for density (kPa)
            std::cout <<"================================================================================";std::cout << "\n";
            std::cout <<"\033[93mAGA 8 GROSS METHOD\033[0m\n";
            std::cout <<"Pressure [MPa.g]     = "; Pin = get_double("LhGreenF"); std::cout << "\n";;
            std::cout <<"Temperature [C]      = "; Tin = get_double("LhGreenF"); std::cout << "\n";;
            std::cout <<"Tref Heating[C]      = "; Th = get_double("LhGreenF"); std::cout << "\n";
            std::cout <<"Tref Density[C]      = "; Td = get_double("LhGreenF"); std::cout << "\n";
            std::cout <<"Pref Density[MPa.g]  = "; Pd = get_double("LhGreenF"); std::cout << "\n";
            TT=Tin+273.15;
            PP=(Pin+atm)*1000.0; //PRESSURE is ABSOLUTE
            Th+=273.15;
            Td+=273.15;
            Pd = (Pd+atm)*1000.0;
            AGA8Gr.MolarMass(mix[activeMix].mix);
            double MW = AGA8Gr.GetMW();
            double Mwl;
            AGA8Gr.GrossHv(HN,HCH);
            printf("HN            = %7.4f [kJ/mol]   HCH       = %7.4f [kJ/mol]\n",HN,HCH);
            AGA8Gr.DensityGross(TT,PP,HCH,D,ierr,herr);
            printf("Density       = %10.7f [kg/m3]\n",D*MW);
            AGA8Gr.PressureGross(TT,D,HCH,PP,Z,ierr,herr);
            printf("Z             =  %9.7f []          Density   = %10.7f [ kg/m3]     Molar Mass = %10.7f [g/mol]\n",Z,D*MW,MW);
            printf("Pressure      = %12.5f [Kpa.a]   HN        = %7.4f    [kJ/mol]    HCH        = %7.4f  [kJ/mol]\n",PP,HN,HCH);
            AGA8Gr.GrossInputs(TT, PP, Gr, HN, HCH, ierr, herr);
            printf("Rel. density  =  %9.7f []          HN        = %7.4f [kJ/mol]    HCH        = %7.4f [kJ/mol]\n",Gr,HN,HCH);
            AGA8Gr.GrossMethod1(Th, Td, Pd, Gr, HN*D, Mwl, HCH, HN, ierr, herr);
            printf("HN  (MET 1)   = %7.4f [kJ/mol]   HCH       = %7.4f [kJ/mol]\n",HN,HCH);
            AGA8Gr.GrossMethod2(Th, Td, Pd, Gr, Hv, Mwl, HCH, HN, ierr, herr);
            printf("HN  (MET 2)   = %7.4f [kJ/mol]   HCH       = %7.4f [kJ/mol]    Hv         = %7.4f [MJ/m^3]\n\n",HN,HCH, Hv);
            choice("Exit (y) ","yY");
        }

        if(charIn=='4' && MixValid){
            std::cout <<"================================================================================";std::cout << "\n";
            std::cout <<"\033[93mGERG 2008\033[0m\n";
            std::cout <<"Pressure [MPa.g] = "; Pin = get_double("LhGreenF"); std::cout << "\n";
            std::cout <<"Temperature [C]  = "; Tin = get_double("LhGreenF"); std::cout << "\n";
            //std::cout <<"Density Flag     = "; dFlag = get_int("LhGreenF"); std::cout << "\n";
            dFlag = static_cast<int> (choice("Density Flag (0/1/2) = ","012"))-48;
            logEOSin(logs, "GERG 2008 METHOD", Pin, Tin, atm);
            TT=Tin+273.15;
            PP=(Pin+atm)*1000.0; //PRESSURE is ABSOLUTE
            Gerg.MolarMass(mix[activeMix].mix);
            double MW = Gerg.GetMW();
            Gerg.Density(dFlag, TT, PP, D, ierr, herr);
            std::cout <<"\033[93mCheck\033[0m\n";
            logEOSdata(logs,"Density Error ID... {:<d}\n", ierr);
            logEOSdata(logs,"Density Error...... {:<s}\n", herr);
            Gerg.Properties(TT, D, P, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa, A);
            std::cout <<"\033[93mResults\033[0m\n";
            if(dFlag==0)logEOSdata(logs,"Density Flag....... {:<s}\n", "strict pressure solver in the gas phase without checks");
            if(dFlag==0)logEOSdata(logs,"                    {:<s}\n", "(fastest mode, but output state may not be stable single phase)");
            if(dFlag==1)logEOSdata(logs,"Density Flag....... {:<s}\n", "checks for possible 2-phase states");
            if(dFlag==1)logEOSdata(logs,"                    {:<s}\n", "(result may still not be stable single phase, but many unstable states will be identified)");
            if(dFlag==2)logEOSdata(logs,"Density Flag....... {:<s}\n", "search for liquid phase (and make the same checks when iFlag=1)");
            logEOSdata(logs,"Molar Mass......... {:<15.4f} {}\n", MW,"[g/mol]");
            logEOSdata(logs,"Density............ {:<15.3f} {}\n", D*MW ,"[kg/m3]"); // [mol/l]*[g/mol]=[g/l]=[Kg/m3]
            logEOSdata(logs,"Pressure........... {:<15.4f} {:10s}       :  {:<15.4f} {}\n", PP/1000.0,"[Mpa.a]", PP*kpaTOpsi, "[psi.a]" );
            logEOSdata(logs,"Z.................. {:<15.6f} {}\n", Z ,"[-]");
            logEOSdata(logs,"U.................. {:<15.2f} {}\n", U/MW, "[KJ/Kg]");
            logEOSdata(logs,"H.................. {:<15.2f} {}\n", H/MW, "[KJ/Kg]");
            logEOSdata(logs,"S.................. {:<15.4f} {}\n", S/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cv................. {:<15.4f} {}\n", Cv/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cp................. {:<15.4f} {}\n", Cp/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Sound Speed........ {:<15.2f} {}\n", W, "[m/S]");
            logEOSdata(logs,"G.................. {:<15.2f} {}\n", G/MW, "[KJ/Kg]");
            logEOSdata(logs,"Joule-Thomson...... {:<15.4f} {}\n", JT*1000.0, "[K/MPa]");
            logEOSdata(logs,"Isentropic Exponent {:<15.7f} {}\n", Kappa, "[-]");
            logEOSdata(logs,"Helmholtz.......... {:<15.4f} {}\n\n", A/MW, "[KJ/Kg]");
            choice("Exit (y) ","yY");
        }

        if(charIn=='5' && MixValid){
            double a0, Da0Dt,Da0Dtt;
            std::cout <<"================================================================================";std::cout << "\n";
            std::cout <<"\033[93mPENG ROBINSON\033[0m\n";
            std::cout <<"Pressure [MPa.g] = "; Pin = get_double("LhGreenF"); std::cout << "\n";
            std::cout <<"Temperature [C]  = "; Tin = get_double("LhGreenF"); std::cout << "\n";
            logEOSin(logs, "PENG ROBINSON METHOD", Pin, Tin, atm);
            TT=Tin+273.15;
            PP=(Pin+atm)*1000.0; //PRESSURE is ABSOLUTE
            Gerg.MolarMass(mix[activeMix].mix);
            Gerg.Density(1, TT, PP, D, ierr, herr);
            Gerg.PropertiesIdeal(TT, D, a0, Da0Dt,Da0Dtt);

            PR.SetIdeal(a0,Da0Dt,Da0Dtt);
            PR.MolarMass(mix[activeMix].mix);
            double MW = PR.GetMW();
            PR.Density(1, TT, PP, D, ierr, herr);
            std::cout <<"\033[93mCheck\033[0m\n";
            logEOSdata(logs,"Density Error ID... {:<d}\n", ierr);
            logEOSdata(logs,"Density Error...... {:<s}\n", herr);
            PR.Properties(TT, D, P, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa, A);
            std::cout <<"\033[93mResults\033[0m\n";
            logEOSdata(logs,"Molar Mass......... {:<15.4f} {}\n", MW,"[g/mol]");
            logEOSdata(logs,"Density............ {:<15.3f} {}\n", D*MW ,"[kg/m3]");
            logEOSdata(logs,"Pressure........... {:<15.4f} {:10s}       :  {:<15.4f} {}\n", PP/1000.0,"[Mpa.a]", PP*kpaTOpsi, "[psi.a]" );
            logEOSdata(logs,"Z.................. {:<15.6f} {}\n", Z ,"[-]");
            logEOSdata(logs,"U.................. {:<15.2f} {}\n", U/MW, "[KJ/Kg]");
            logEOSdata(logs,"H.................. {:<15.2f} {}\n", H/MW, "[KJ/Kg]");
            logEOSdata(logs,"S.................. {:<15.4f} {}\n", S/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cv................. {:<15.4f} {}\n", Cv/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cp................. {:<15.4f} {}\n", Cp/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Sound Speed........ {:<15.2f} {}\n", W, "[m/S]");
            logEOSdata(logs,"G.................. {:<15.2f} {}\n", G/MW, "[KJ/Kg]");
            logEOSdata(logs,"Joule-Thomson...... {:<15.4f} {}\n", JT*1000.0, "[K/MPa]");
            logEOSdata(logs,"Isentropic Exponent {:<15.7f} {}\n", Kappa, "[-]");
            logEOSdata(logs,"Helmholtz.......... {:<15.4f} {}\n\n", A/MW, "[KJ/kg]");
            choice("Exit (y) ","yY");
        }

        if(charIn=='6' && MixValid){
            std::cout <<"================================================================================";std::cout << "\n";
            std::cout <<"\033[93mBWRS\033[0m\n";
            std::cout <<"Pressure [MPa.g] = "; Pin = get_double("LhGreenF"); std::cout << "\n";
            std::cout <<"Temperature [C]  = "; Tin = get_double("LhGreenF"); std::cout << "\n";
            logEOSin(logs, "BWRS", Pin, Tin, atm);
            TT=Tin+273.15;
            PP=(Pin+atm); //MPa PRESSURE is ABSOLUTE
            Bwrs.MolarMass(mix[activeMix].mix);
            double MW = Bwrs.GetMW();
            Bwrs.Density(0, TT, PP, D, ierr, herr); // density g/m3
            std::cout <<"\033[93mCheck\033[0m\n";
            logEOSdata(logs,"Density Error ID... {:<d}\n", ierr);
            logEOSdata(logs,"Density Error...... {:<s}\n", herr);
            Bwrs.Properties(TT, D, P, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa, A);
            std::cout <<"\033[93mResults\033[0m\n";
            logEOSdata(logs,"Molar Mass......... {:<15.4f} {}\n", MW,"[g/mol]");
            logEOSdata(logs,"Density............ {:<15.3f} {}\n", D*1000.0*MW ,"[kg/m3]"); // convert density from from g/m3 tp Kg/m3
            logEOSdata(logs,"Pressure........... {:<15.4f} {:10s}       :  {:<15.4f} {}\n", PP,"[Mpa.a]", PP*kpaTOpsi, "[psi.a]" );
            logEOSdata(logs,"Z.................. {:<15.6f} {}\n", Z ,"[-]");
            logEOSdata(logs,"U.................. {:<15.2f} {}\n", U/MW, "[KJ/Kg]");
            logEOSdata(logs,"H.................. {:<15.2f} {}\n", H/MW, "[KJ/Kg]");
            logEOSdata(logs,"S.................. {:<15.4f} {}\n", S/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cv................. {:<15.4f} {}\n", Cv/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cp................. {:<15.4f} {}\n", Cp/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Sound Speed........ {:<15.2f} {}\n", W, "[m/S]");
            logEOSdata(logs,"G.................. {:<15.2f} {}\n", G/MW, "[KJ/Kg]");
            logEOSdata(logs,"Joule-Thomson...... {:<15.4f} {}\n", JT*1000.0, "[K/MPa]");
            logEOSdata(logs,"Isentropic Exponent {:<15.7f} {}\n", Kappa, "[-]");
            logEOSdata(logs,"Helmholtz.......... {:<15.4f} {}\n\n", A/MW, "[KJ/kg]");
            choice("Exit (y) ","yY");
        }

        if(charIn=='9'){
            std::cout <<"================================================================================";std::cout << "\n";
            std::cout <<"\033[93mIAPWS 97 [2012]\033[0m\n";
            std::cout <<"Pressure [MPa.g] = "; Pin = get_double("LhGreenF"); std::cout << "\n";
            std::cout <<"Temperature [C]  = "; Tin = get_double("LhGreenF"); std::cout << "\n";
            logEOSin(logs, "IAPWS 97[2012]", Pin, Tin, atm);
            TT=Tin+273.15;
            PP=Pin+atm; //PRESSURE is ABSOLUTE
            std::cout <<"\033[93mResults\033[0m\n";
            logEOSdata(logs,"Rho................ {:<15.4f} {}\n", IF97::rhomass_Tp(TT, PP),"[Kg/m3]");
            logEOSdata(logs,"Rho................ {:<15.11f} {}\n", 1.0/IF97::rhomass_Tp(TT, PP),"[m3/Kg]");
            logEOSdata(logs,"U.................. {:<15.5f} {}\n", IF97::umass_Tp(TT, PP), "[KJ/Kg]");
            logEOSdata(logs,"H.................. {:<15.5f} {}\n", IF97::hmass_Tp(TT, PP), "[KJ/Kg]");
            logEOSdata(logs,"S.................. {:<15.8f} {}\n", IF97::smass_Tp(TT, PP), "[KJ/Kg K]");
            logEOSdata(logs,"Cv................. {:<15.7f} {}\n", IF97::cvmass_Tp(TT, PP), "[KJ/Kg K]");
            logEOSdata(logs,"Cp................. {:<15.7f} {}\n", IF97::cpmass_Tp(TT, PP), "[KJ/Kg K]");
            logEOSdata(logs,"w.................. {:<15.6f} {}\n", IF97::speed_sound_Tp(TT, PP), "[m/s]");
            choice("Exit (y) ","yY");
        }

    if(charIn=='r' && MixValid){
        std::cout <<"================================================================================";std::cout << "\n";
        double ts, te,td, ps, pe, pd;
        std::cout <<  "Temperature START    [\370C] "; ts=get_double("LhGreenF");
        std::cout <<"\nTemperature END      [\370C] "; te=get_double("LhGreenF");
        std::cout <<"\nTemperature STEP     [\370C] "; td=get_double("LhGreenF");
        std::cout <<"\n\nPressure    START [Mpa.g] "; ps=get_double("LhGreenF");
        std::cout <<  "\nPressure    END   [Mpa.g] "; pe=get_double("LhGreenF");
        std::cout <<  "\nPressure    STEP  [Mpa.g] "; pd=get_double("LhGreenF");

        std::cout <<"\n\nFile name without extension "; fname=get_string("LhGreenF");
        if(fname.empty()) fname="test";
        std::cout <<'\n';
        std::filesystem::path directorypath;
        directorypath = std::filesystem::current_path();
        //fname = directorypath.string().c_str()+fname+".csv";
        fname = fname+".csv";
        FILE * out;
        if((out = fopen(fname.c_str(), "w")) == NULL){
            std::cout << "Could not open the file - '" << fname << "'\n";
        } else {
            printf("\n            T[\370C] ");
            fprintf(out,";T[\370C]; ");
            for (double ttemp = ts; ttemp<=te; ttemp+=td){ //PRINT TEMPERATURE
                printf("%6.2f     ",ttemp);
                fprintf(out,"%6.2f; ",ttemp); // write
            }
            printf("\n");
            fprintf(out,"\n"); // write
            AGA8Det.MolarMass(mix[activeMix].mix);
            MW = AGA8Det.GetMW();

            for (; ps<=pe; ps+=pd){ //PRESSSUE
                printf("P[Mpa.g] %7.2f ",ps);
                fprintf(out,"P[Mpa.g]; %7.2f; ",ps);
                for (double ttemp = ts; ttemp<=te; ttemp+=td){ //TEMPERATURE
                    //for (size_t k=0; k<mix.size(); k++){ //if you want to cycle on mix
                        TT=ttemp+273.15;
                        PP=ps+atm; //PRESSURE is ABSOLUTE

                        D = -1.0;
                        AGA8Det.Density(TT, PP*1000.0 ,D, ierr, herr);
                        if(ierr==0) {
                            AGA8Det.Properties(TT, D, PP, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa,A);
                            //fprintf(out,"%9.5f; ",JT*1000.0);
                            //printf("%9.5f; ",JT*1000.0);
                            fprintf(out,"%9.5f; ",D*MW);
                            printf("%9.5f; ",D*MW);
                            //fprintf(out,"%15.7f; ",Kappa);
                            //printf("%15.7f; ",Kappa);
                        }

                        //Aga8.DZOFPT(&PP, &TT, &D, &Z, &BMIX);
                        //printf("%9.5f ",Z);
                        //fprintf(out,"%9.5f; ",Z); // write
                        //};
                    };
                printf("\n");
                fprintf(out,"\n"); // write
            };
            fclose(out); // close file
            std::cout << "Task Finish - '" << mix.size() << "'\n";
        }
        choice("Exit (y) ","yY");
    }

    } while (charIn != 'e');
}

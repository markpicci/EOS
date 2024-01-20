#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <stdlib.h>
#include <conio.h>
#include <filesystem>
#include "conIO.h"
#include "ISO-20765-5.hpp"
#include "ISO-12213-2.hpp"
#include "AGA8 Gross.hpp"
#include "AGA8 Detail.hpp"
#include "GERG2008.hpp"
#include "generalizedcubic.hpp"
#include "IF97.hpp"

#define naga 21

std::vector<std::string>  gasName = {"Methane", "Nitrogen", "Carbon_dioxide", "Ethane", "Propane", "iso-Butane", "n-Butane", "iso-Pentane", "n-Pentane", "n-Hexane",
                     "n-Heptane", "n-Octane", "n-Nonane", "n-Decane", "Hydrogen", "Oxygen", "Carbon_monoxide",  "Water", "Hydrogen_sulfide", "Helium", "Argon"};
struct mixSt{
    std::string name;
    std::map<std::string, double> mix;
};
std::ofstream logs;
std::vector<mixSt> mix;
size_t activeMix = -1;
bool MixValid(){ return (activeMix >= 0 &&  activeMix < mix.size());}

template<typename... Args>
void logEOSdata(std::ofstream& logs, std::string_view fmt_str, Args&&... args){
    std::string outFor = std::vformat(fmt_str, std::make_format_args(args...));
    std::cout << outFor;
    logs << outFor;
}

int manageMix(){
std::string fname;
mixSt tmix;
char charIn;
double perc, tot;

    system("cls");
    do{
        std::vector<std::string> menuList = {"(g)as","(l)oad","(s)ave","(i)nput","(m)ix","(c)ollection","(d)elete", "(e)xit"};
        charIn = conMenu("GAS MISTURE",menuList);
        std::cout <<"\n";

        if(charIn=='g'){
            std::cout <<"================================================================================";std::cout << "\n";
            std::filesystem::path directorypath;
            directorypath = std::filesystem::current_path();
            for (const auto & entry : std::filesystem::directory_iterator(directorypath))
                if(entry.path().extension() == ".mix")
                    std::cout << entry.path() << "\n";
            choice("\nReturn Main Menu (y) ","yY");
        }

        if(charIn=='l'){
            std::cout <<"================================================================================";std::cout << "\n";
            std::cout <<"File name without extension "; fname=get_string("LhGreenF");
            if(fname.empty()) {fname="test"; std::cout <<fname;}
            std::cout <<'\n';
            std::filesystem::path directorypath;
            directorypath = std::filesystem::current_path();
            std::ifstream fh(fname+".mix");
            if( fh.is_open()){
                int nmix;
                fh >> nmix;
                for(int i=0;i<nmix;i++){
                    fh>>tmix.name;
                    for(int i =0; i<naga;i++){
                        std::string gasName; double gasComp;
                        fh>>gasName;
                        fh>>gasComp;
                        tmix.mix[gasName] = gasComp;
                    }
                    mix.push_back(tmix); activeMix = mix.size()-1;
                }
                fh.close(); // close file
                std::cout << "Number of mixture loaded - '" << mix.size() << "'\n";

                if (activeMix==0){
                    logMix(logs,mix[activeMix].name, mix[activeMix].mix);
                    choice("Exit (y) ","yY");
                } else {
                    int i=0;
                    std::cout <<"\nSelect Active Mixture\n";
                    for(const auto& mixPair:mix) {printf("%2d     Mixture Name     \033[92m%s\033[0m\n", i, mixPair.name.c_str());i++;}
                    std::cout <<"\nSelect Active Mixture "; activeMix = get_int("LhYellowF");std::cout <<"\n";
                    if(!MixValid()){
                        activeMix = -1;
                        std::cout << "Selection NON valid \n";
                        choice("Exit (y) ","yY");
                    } else
                        logMix(logs,mix[activeMix].name, mix[activeMix].mix);
                }
            } else {std::cout << "Could not open the file - '" << fname << "'\n"; choice("Exit (y) ","yY");}
        }

        if(charIn=='s'){
            std::cout <<"================================================================================";std::cout << "\n";
            std::cout <<"File name without extension "; fname=get_string("LhGreenF");
            if(fname.empty()) {fname="test"; std::cout <<fname;}
            std::cout <<'\n';
            std::filesystem::path directorypath;
            directorypath = std::filesystem::current_path();
            std::ofstream fh(fname+".mix");
            if( fh.is_open()){
                int nmix=mix.size();
                fh << nmix<<'\n';
                for(int i=0;i<nmix;i++){
                    fh<<mix[i].name<<'\n';
                    for(const auto& gasPair: mix[i].mix){
                        fh<<gasPair.first<<" ";
                        fh<<gasPair.second<<'\n';
                    }
                }
                fh.close(); // close file
                std::cout << "Number of mixture saved - '" << mix.size() << "'\n";
                getch();
            } else std::cout << "Could not open the file - '" << fname << "'\n";
        }

        if(charIn=='d'){
            std::cout <<"================================================================================";std::cout << "\n";
            charIn=choice("Confirm to DELETE all mixtures (y/n) ","ynYN");
            if (charIn =='y' || charIn =='Y'){
                    mix.clear();
                    activeMix = -1;
            }
        }

        if(charIn=='i'){
            system("cls");
            std::cout <<"================================================================================";std::cout << "\n";
            std::cout <<"Mixture Name  "; tmix.name=get_string(); std::cout << '\n';
            tot = 0.0;
            for(int i =0; i<naga;i++){
                printf("%2d  \033[92m%18s\033[0m ",i, gasName[i].c_str());
                perc=get_double(); tmix.mix[gasName[i]] = perc;
                printf("\n");
                tot += perc;
            }
            printf("\n");
            charIn=choice("Confirm to add (y/n) ","ynYN");
            printf("\n       Total composition %7.5f\n\n", tot);
            if (charIn =='y' || charIn =='Y') {
                mix.push_back(tmix);
                activeMix = mix.size()-1;
                logMix(logs,mix[activeMix].name, mix[activeMix].mix);
            }
        }

        if(charIn=='m' && !MixValid()){
            std::cout <<"================================================================================";std::cout << "\n";
            choice("NOT valid mixture       Exit (y) \n","yY");printf("\n");
        }
        if(charIn=='m' && MixValid()){
            system("cls");
            std::cout <<"================================================================================";std::cout << "\n";
            tmix = mix[activeMix];
            printf("      Mixture Name     \033[93m%s\033[0m\n", tmix.name.c_str());
            tot = 0.0;
            for(int i =0; i<naga;i++){
                tot += tmix.mix[gasName[i]];
                if (tmix.mix[gasName[i]]!= 0.0)
                    printf("%2d  \033[92m%18s\033[0m xi= %7.5f\n",i, gasName[i].c_str(), tmix.mix[gasName[i]]);
            }
            printf("\n       Total composition %7.5f\n", tot);
            choice("       Exit (y) ","yY");
        }

        if(charIn=='c'){
            int i=0;
            system("cls");
            std::cout <<"================================================================================";std::cout << "\n";
            for(const auto& mixPair:mix) {printf("%2d     Mixture Name     \033[92m%s\033[0m\n", i, mixPair.name.c_str());i++;}
            std::cout <<"\nSelect Active Mixture "; activeMix = get_int("LhYellowF");
            if(!MixValid())
                activeMix = -1;
            else
                logMix(logs,mix[activeMix].name, mix[activeMix].mix);
        }

    } while (charIn != 'e');
    return 0;
}

int main(int argc, char *argv[]) {

double TT, PP, BMIX;
double D, P, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa, A;

int ierr, dFlag;
std::string herr;
std::string fname;
FILE* out;
ISO122132 Aga8;
ISO207655 Visosity;
AGA8Gross  AGA8Gr;
AGA8Detail AGA8Det;
GERG2008   Gerg;
PengRobinson PR;


double Tin, Pin, MW;
char charIn;

    initCon();
    std::cout <<"LOG File name without extension "; fname=get_string("LhGreenF");
    if(fname.empty()) {fname="log"; std::cout <<fname;}
    std::cout <<'\n';
    std::filesystem::path directorypath;
    directorypath = std::filesystem::current_path();
    logs.open(fname+".log");
    if (!logs.is_open()){
        std::cout <<"ERROR opening LOG file\n";
        return 1;
    }
    logHeader(logs);
    std::vector<std::string> menuList = {"(1) ISO-12213-2","(2) AGA8 Detail","(3) AGA8Gross","(4) GERG 2008","(5) PENG ROBINSON","(9) IAPWS R7-97 [2012]", "", "(m)ix","(r)ange","(e)xit"};

//===================================== MANUAL
    do{
        charIn = conMenu("MAIN MENU",menuList);
        std::cout <<"\n";

        if(charIn=='m') manageMix();

        if(charIn=='1' && MixValid()){
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
            logEOSdata(logs,"Density............ {:<15.2f} {}\n", D*MW,"[kg/m3]");
            logEOSdata(logs,"Z.................. {:<15.6f} {}\n", Z   ,"[-]");
            logEOSdata(logs,"Viscosity.......... {:<15.4f} {}\n\n", Visosity.VIS(TT, D, mix[activeMix].mix),"[mPa*s]");
            choice("Exit (y) ","yY");
        }

        if(charIn=='2' && MixValid()){
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
            logEOSdata(logs,"Density............ {:<15.2f} {}\n", D*MW ,"[kg/m3]");
            logEOSdata(logs,"Pressure........... {:<15.4f} {:10s}       :  {:<15.4f} {}\n", PP/1000.0,"[Mpa.a]", PP*kpaTOpsi, "[psi.a]" );
            logEOSdata(logs,"Z.................. {:<15.6f} {}\n", Z ,"[-]");
            logEOSdata(logs,"U.................. {:<15.2f} {}\n", U/MW, "[KJ/Kg]");
            logEOSdata(logs,"H.................. {:<15.2f} {}\n", H/MW, "[KJ/Kg]");
            logEOSdata(logs,"S.................. {:<15.4f} {}\n", S/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cv................. {:<15.4f} {}\n", Cv/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cp................. {:<15.4f} {}\n", Cp/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Sound Speed........ {:<15.2f} {}\n", W, "[m/S]");
            logEOSdata(logs,"G.................. {:<15.2f} {}\n", G/MW, "[KJ/Kg]");
            logEOSdata(logs,"Joule-Thomson...... {:<15.5f} {}\n", JT, "[K/kPa]");
            logEOSdata(logs,"Isentropic Exponent {:<15.7f} {}\n", Kappa, "[-]");
            logEOSdata(logs,"Helmholtz.......... {:<15.4f} {}\n\n", A/MW, "[KJ/kPa]");
            choice("Exit (y) ","yY");
        }

        if(charIn=='3' && MixValid()){
            double Hv, Gr, HN, HCH;
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

        if(charIn=='4' && MixValid()){
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
            logEOSdata(logs,"Density Flag....... {:<d}\n", dFlag);
            logEOSdata(logs,"Density Error ID... {:<d}\n", ierr);
            logEOSdata(logs,"Density Error...... {:<s}\n", herr);
            Gerg.Properties(TT, D, P, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa, A);
            std::cout <<"\033[93mResults\033[0m\n";
            logEOSdata(logs,"Molar Mass......... {:<15.4f} {}\n", MW,"[g/mol]");
            logEOSdata(logs,"Density............ {:<15.2f} {}\n", D*MW ,"[kg/m3]"); // [mol/l]*[g/mol]=[g/l]=[Kg/m3]
            logEOSdata(logs,"Pressure........... {:<15.4f} {:10s}       :  {:<15.4f} {}\n", PP/1000.0,"[Mpa.a]", PP*kpaTOpsi, "[psi.a]" );
            logEOSdata(logs,"Z.................. {:<15.6f} {}\n", Z ,"[-]");
            logEOSdata(logs,"U.................. {:<15.2f} {}\n", U/MW, "[KJ/Kg]");
            logEOSdata(logs,"H.................. {:<15.2f} {}\n", H/MW, "[KJ/Kg]");
            logEOSdata(logs,"S.................. {:<15.4f} {}\n", S/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cv................. {:<15.4f} {}\n", Cv/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cp................. {:<15.4f} {}\n", Cp/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Sound Speed........ {:<15.2f} {}\n", W, "[m/S]");
            logEOSdata(logs,"G.................. {:<15.2f} {}\n", G/MW, "[KJ/Kg]");
            logEOSdata(logs,"Joule-Thomson...... {:<15.7f} {}\n", JT, "[K/kPa]");
            logEOSdata(logs,"Isentropic Exponent {:<15.7f} {}\n", Kappa, "[-]");
            logEOSdata(logs,"Helmholtz.......... {:<15.4f} {}\n\n", A/MW, "[KJ/Kg]");
            choice("Exit (y) ","yY");
        }

        if(charIn=='5' && MixValid()){
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
            logEOSdata(logs,"Density............ {:<15.2f} {}\n", D*MW ,"[kg/m3]");
            logEOSdata(logs,"Pressure........... {:<15.4f} {:10s}       :  {:<15.4f} {}\n", PP/1000.0,"[Mpa.a]", PP*kpaTOpsi, "[psi.a]" );
            logEOSdata(logs,"Z.................. {:<15.6f} {}\n", Z ,"[-]");
            logEOSdata(logs,"U.................. {:<15.2f} {}\n", U/MW, "[KJ/Kg]");
            logEOSdata(logs,"H.................. {:<15.2f} {}\n", H/MW, "[KJ/Kg]");
            logEOSdata(logs,"S.................. {:<15.4f} {}\n", S/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cv................. {:<15.4f} {}\n", Cv/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Cp................. {:<15.4f} {}\n", Cp/MW, "[KJ/(Kg-K)]");
            logEOSdata(logs,"Sound Speed........ {:<15.2f} {}\n", W, "[m/S]");
            logEOSdata(logs,"G.................. {:<15.2f} {}\n", G/MW, "[KJ/Kg]");
            logEOSdata(logs,"Joule-Thomson...... {:<15.7f} {}\n", JT, "[K/kPa]");
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

    if(charIn=='r' && MixValid()){
        std::cout <<"================================================================================";std::cout << "\n";
        double ts, te,td, ps, pe, pd;
        std::cout <<  "Temperature START    [\370C] "; ts=get_double("LhGreenF");
        std::cout <<"\nTemperature END      [\370C] "; te=get_double("LhGreenF");
        std::cout <<"\nTemperature STEP     [\370C] "; td=get_double("LhGreenF");
        std::cout <<"\n\nPressure    START [Mpa.g] "; ps=get_double("LhGreenF");
        std::cout <<  "\nPressure    END   [Mpa.g] "; pe=get_double("LhGreenF");
        std::cout <<  "\nPressure    STEP  [Mpa.g] "; pd=get_double("LhGreenF");

        std::cout <<"\n\nFile name without extension "; fname=get_string("LhGreenF");
        if(fname.empty()) {fname="test"; std::cout <<fname;}
        std::cout <<'\n';
        std::filesystem::path directorypath;
        directorypath = std::filesystem::current_path();
        fname = directorypath.string().c_str()+fname+".csv";
        if((out = fopen(fname.c_str(), "w")) == NULL){
            std::cout << "Could not open the file - '" << fname << "'\n";
        } else {
            printf("\n            T[\370C] ");
            fprintf(out,"T[\370C]; ");
            for (double ttemp = ts; ttemp<=te; ttemp+=td){ //PRINT TEMPERATURE
                printf("%6.2f     ",ttemp);
                fprintf(out,"%6.2f; ",ttemp); // write
            }
            printf("\n");
            fprintf(out,"\n"); // write
            Aga8.MolarMass(mix[activeMix].mix);

            for (; ps<=pe; ps+=pd){ //PRESSSUE
                printf("P[Mpa.g] %7.2f ",ps);
                fprintf(out,"P[Mpa.g]; %7.2f; ",ps);
                for (double ttemp = ts; ttemp<=te; ttemp+=td){ //TEMPERATURE
                    { for (size_t k=0; k<mix.size(); k++){
                        TT=ttemp+273.15;
                        PP=ps+atm; //PRESSURE is ABSOLUTE
                        Aga8.DZOFPT(&PP, &TT, &D, &Z, &BMIX);
                        printf("%9.5f ",Z);
                        fprintf(out,"%9.5f; ",Z); // write
                        };
                    };
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
    logs.close();
}

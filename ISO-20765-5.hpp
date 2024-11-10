#ifndef __ISO20765_5__H
#define __ISO20765_5__H

#include <math.h>
//  **********  This code is preliminary, and will be updated.  **********
//  **********  Use only for beta testing.  **********

// ISO 12213-2: Pipeline gas compressibility computation according to the method AGA8-92DC
class ISO207655{
private:
//===============================================================================================================
//ISO-20765-5-2022
//dynamic viscosity mPa·s
//
// Temperature K
// Density Molar density mol/dm3
//
enum {naga = 21};
//std::string gas[naga] = {"Methane",          "Nitrogen", "Carbon dioxide",  "Ethane",    "Propane",  "Water",
                         //"Hydrogen sulfide", "Hydrogen", "Carbon monoxide", "Oxygen",
                         //"iso-Butane",       "n-Butane", "iso-Pentane",     "n-Pentane", "n-Hexane",
                         //"n-Heptane",        "n-Octane", "n-Nonane",        "n-Decane",  "Helium",   "Argon"};

std::vector<std::string> gas = {"74-82-8" ,        "7727-37-9",       "124-38-9",         "74-84-0",      "74-98-6",   "7732-18-5",
                                "7783-06-4",        "1333-74-0",       "630-08-0",         "7782-44-7",
                                "75-28-5",          "106-97-8",        "78-78-4",          "109-66-0",     "110-54-3",
                                "142-82-5",         "111-65-9",        "111-84-2",         "124-18-5",     "7440-59-7", "7440-37-1"};

double visMi[naga] = {16.04246, 28.0134, 44.0095, 30.06904, 44.09562, 18.01528, 34.08088, 2.01588, 28.0101, 31.9988, 58.1222, 58.1222, 72.14878, 72.14878, 86.17536, 100.20194, 114.22852, 128.2551, 142.28168, 4.002602, 39.948}; // (g*mol-1)
double visRhoc[naga] = {10.139342719, 11.1839, 10.624978698, 6.87085454, 5.000043088, 17.87371609, 10.19, 14.94, 10.85, 13.63, 3.86014294, 3.920016792, 3.271, 3.215577588, 2.705877875, 2.315324434, 2.056404127, 1.81, 1.64, 17.399, 13.407429659}; //(mol*dm-3)
double visTc[naga] = {190.564, 126.192, 304.1282, 305.322, 369.825, 647.096, 373.1, 33.19, 132.86, 154.595, 407.817, 425.125, 460.35, 469.7, 507.82, 540.13, 569.32, 594.55, 617.7, 5.1953, 150.687}; // (K)
double visPc[naga] = {4.5992, 3.3958, 7.3773, 4.8722, 4.2512, 22.064, 9, 1.2964, 3.494, 5.043, 3.629, 3.796, 3.378, 3.37, 3.034, 2.736, 2.497, 2.281, 2.103, 0.2276, 4.863}; // (MPa)

public:

    ISO207655(){}
    ~ISO207655(){}

    double VIS(double &T, double &D, std::map<std::string, double>& mix){
    int   I;
    double visi[naga], XJ[naga];
    double Mmix, Tcmix, Pcmix, Vcmix, uNeta, uP, Tri, rhorMix, alfa, delta;

        for(std::size_t i = 1; i < naga; ++i) XJ[i] = mix[gas[i]];

        Mmix=Tcmix=Pcmix=Vcmix =0.0;
        for(I=0; I<naga; I++){
            Mmix += XJ[I]*visMi[I];
            Tcmix += XJ[I]*visTc[I];
            Pcmix += XJ[I]*visPc[I];
            Vcmix += XJ[I]/visRhoc[I];
        }

        uNeta = 0.0001; //mPa·s
        uP = 0.101325; //MPa

        for(I=0; I<naga; I++){
            visi[I]=uNeta*pow(visMi[I],0.5)*pow(visTc[I],-1.0/6.0)*pow(visPc[I]/uP,2.0/3.0);
            Tri = T / visTc[I];
            if(Tri<=1.5){
                alfa = 3.4*pow(Tri,0.94);
            }else{
                alfa = 1.778*pow(4.58*Tri-1.67,0.625);
            }
            visi[I] *= alfa;
        }

        rhorMix = Vcmix* D; //where D is the molar density at T and P, calculated from ISO 20765-2.
        delta = 1.023+0.23364*rhorMix+0.58533*pow(rhorMix,2)-0.40758*pow(rhorMix,3)+0.093324*pow(rhorMix,4);

    double num, den, visMix,visFin, eta;
        num=den=0.0;
        for(I=0; I<naga; I++){
            num += XJ[I]*visi[I]*pow(visMi[I],0.5);
            den += XJ[I]*pow(visMi[I],0.5);
        }
        visMix = num/den;
        eta = uNeta*pow(Mmix,0.5)*pow(Tcmix,-1.0/6.0)*pow(Pcmix/uP,2.0/3.0);
        visFin = visMix+eta*(pow(delta,4)-1.0);
        return visFin;
    }
};

#endif // ISO20765

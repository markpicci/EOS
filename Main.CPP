#include "mainh.hpp"
#include "conIO.h"
#include "DB.hpp"

#define naga 21

DB db;
//std::vector<std::string>  gasName = {"Methane",    "Nitrogen",        "Carbon_dioxide", "Ethane",           "Propane",
                                     //"iso-Butane", "n-Butane",        "iso-Pentane",    "n-Pentane",        "n-Hexane",
//                                     "n-Heptane",  "n-Octane",        "n-Nonane",       "n-Decane",         "Hydrogen",
                                    // "Oxygen",     "Carbon_monoxide", "Water",          "Hydrogen_sulfide", "Helium",   "Argon"};

std::vector<std::string>  gasName =   {"74-82-8",    "7727-37-9",       "124-38-9",       "74-84-0",          "74-98-6",
                                       "75-28-5",    "106-97-8",        "78-78-4",        "109-66-0",         "110-54-3",
                                       "142-82-5",   "111-65-9",        "111-84-2",       "124-18-5",         "1333-74-0",
                                       "7782-44-7",  "630-08-0",        "7732-18-5",      "7783-06-4",        "7440-59-7", "7440-37-1"};


template<typename... Args>
void logEOSdata(std::ofstream& logs, std::string_view fmt_str, Args&&... args){
    std::string outFor = std::vformat(fmt_str, std::make_format_args(args...));
    std::cout << outFor;
    logs << outFor;
}

int manageMix(std::ofstream& logs, size_t& activeMix, std::vector<mixSt>& mix){
std::string fname;
mixSt tmix;
char charIn;
double perc, tot;
bool MixValid;

    system("cls");
    do{
        std::vector<std::string> menuList = {"(g)as","(l)oad","(s)ave","(i)nput","(m)ix","(c)ollection","(d)elete", "(e)xit"};
        charIn = conMenu("GAS MISTURE",menuList);
        std::cout <<"\n";
        MixValid = (activeMix >= 0 &&  activeMix < mix.size());

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
                    MixValid = (activeMix >= 0 &&  activeMix < mix.size());
                    if(!MixValid){
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
            do{
                system("cls");
                std::cout <<"================================================================================";std::cout << "\n";
                std::cout <<"Mixture Name  "; tmix.name=get_string(); std::cout << "\n\n";
                tot = 0.0;
                for(int i =0; i<naga;i++){
                    bool res;
                    std::string gasN ;
                    gasN = std::any_cast<std::string>(db.get(gasName[i], DB::Eq::GEN, DB::DataID::name, res));
                    printf("%2d  \033[92m%18s\033[0m ",i, gasN.c_str());
                    perc=get_double(); tmix.mix[gasName[i]] = perc;
                    printf("\n");
                    tot += perc;
                }
                printf("\n");
                charIn=choice("Confirm (y) or Revise (r)","yYrR");
            } while(charIn !='y' && charIn !='Y');
            printf("\n");
            charIn=choice("Confirm to add (y/n) ","ynYN");
            printf("\n       Total composition %7.5f\n\n", tot);
            if (charIn =='y' || charIn =='Y') {
                mix.push_back(tmix);
                activeMix = mix.size()-1;
                logMix(logs,mix[activeMix].name, mix[activeMix].mix);
            }
        }

        if(charIn=='m' && !MixValid){
            std::cout <<"================================================================================";std::cout << "\n";
            choice("NOT valid mixture       Exit (y) \n","yY");printf("\n");
        }
        if(charIn=='m' && MixValid){
            system("cls");
            std::cout <<"================================================================================";std::cout << "\n";
            tmix = mix[activeMix];
            printf("      Mixture Name     \033[93m%s\033[0m\n\n", tmix.name.c_str());
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
            MixValid = (activeMix >= 0 &&  activeMix < mix.size());
            if(!MixValid)
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
double Tin, Pin, MW;
char charIn;
int ierr, dFlag;
DB db;
std::string herr;
std::string fname;
std::ofstream logs;
std::vector<mixSt> mix;
size_t activeMix = -1;

    MaximizeWindow();
    clearScreen();
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
    std::vector<std::string> menuList = {"(1) EOS", "(2) equipment", "(3) mix","(e)xit"};

    do{
        charIn = conMenu("MAIN MENU",menuList);
        std::cout <<"\n";
        if(charIn=='1') EOS(logs, activeMix, mix);
        if(charIn=='2') equipmentMain(logs);
        if(charIn=='3') manageMix(logs, activeMix, mix);

    } while (charIn != 'e');

    logs.close();
    return 0;
}

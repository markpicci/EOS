#ifndef __SUS_DB_ISA__
#define __SUS_DB_ISA__

#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <any>
#include <cassert>

class DB{
public:
    enum Eq{
        GEN,
        AGA8,
        ISO
    };

    enum DataID{
        CAS,        // CAS number
        name,       // Name of substance
        formula,    // chemical formula
        Mm,         // Mm - Molar mass (g/mol)
        Tc,         // Critical temperatures (in K)
        Pc,         // Critical pressures (in KPa)
        Rhoc,       // Critical Density (in mol/L)
        Zc,         // Critical Compression Factor [-]
        Acentric   // acentric factors (unitless)
    };

private:

    enum {nSos = 21}; //for set up as per AGA of ISO
    struct data{
        double val;
        std::vector<double> vec;
        std::string str;
    };
    std::vector<double> vecEmpty;
    std::string strEmpty = {};
    std::map<std::tuple<std::string, Eq, DataID>, data> DBm;

public:

    ~DB(){}

    DB(){ initSTD();}

    std::any get(std::string CAS_, Eq Eq_, DataID DataID_, bool& findProp){
        std::any ret;
        std::map<std::tuple<std::string, Eq, DataID>, data>::iterator it = DBm.find(std::make_tuple(CAS_, Eq_, DataID_));
        if (it != DBm.end()) {
            findProp = true;
            if(DataID_<3)  ret = (it->second).str;
            if(DataID_>=3) ret = (it->second).val;
            return ret;
        }
        assert(("DB.GET not found", false));
        return 0;
    }

    double getDouble(std::string CAS_, Eq Eq_, DataID DataID_, bool& findProp){
        std::map<std::tuple<std::string, Eq, DataID>, data>::iterator it = DBm.find(std::make_tuple(CAS_, Eq_, DataID_));
        if (it != DBm.end()) {
            findProp = true;
            return (it->second).val ;
        }
        findProp = false;
        return 0.0;
    }

    const std::string& getStr(std::string CAS_, Eq Eq_, DataID DataID_, bool& findProp){
        std::map<std::tuple<std::string, Eq, DataID>, data>::iterator it = DBm.find(std::make_tuple(CAS_, Eq_, DataID_));
        if (it != DBm.end()) {
            findProp = true;
            return (it->second).str;
        }
        findProp = true;
        return strEmpty;
    }

    const std::vector<double>& getVec(std::string CAS_, Eq Eq_, DataID DataID_, bool& findProp){
        std::map<std::tuple<std::string, Eq, DataID>, data>::iterator it = DBm.find(std::make_tuple(CAS_, Eq_, DataID_));
        if (it != DBm.end()) {
            findProp = true;
            return (it->second).vec;
        }
        findProp = 0.0;
        return vecEmpty;
    }


    void initSTD(){
        std::vector<std::string> name_= {"Methane",         "Nitrogen",        "Carbon_dioxide",    "Ethane",       "Propane",  "Water",
                                        "Hydrogen_sulfide", "Hydrogen",        "Carbon_monoxide",   "Oxygen",
                                        "iso-Butane",       "n-Butane",        "iso-Pentane",       "n-Pentane",    "n-Hexane",
                                        "n-Heptane",        "n-Octane",        "n-Nonane",          "n-Decane",     "Helium",   "Argon"};

        std::vector<std::string> CAS_ = {"74-82-8" ,        "7727-37-9",       "124-38-9",         "74-84-0",      "74-98-6",   "7732-18-5",
                                        "7783-06-4",        "1333-74-0",       "630-08-0",         "7782-44-7",
                                        "75-28-5",          "106-97-8",        "78-78-4",          "109-66-0",     "110-54-3",
                                        "142-82-5",         "111-65-9",        "111-84-2",         "124-18-5",     "7440-59-7", "7440-37-1"};

        std::vector<std::string> formula_ = {"CH4" ,        "N2" ,             "CO2" ,              "C2H6",         "C3H8",     "H2O" ,
                                        "H2S" ,             "H2",              "CO",                "O2" ,
                                        "C4H10",            "C4H10",           "C5H12",             "C5H12",        "C6H14" ,
                                        "C7H16",            "C8H18",           "C9H20",             "C10H22",       "He" ,      "Ar" };

        std::vector<std::string> struct_ = {"CH4" ,        "N2" ,             "CO2" ,              "CH3CH3",       "CH3CH2CH3","H2O" ,
                                        "H2S" ,             "H2",              "CO",                "O2" ,
                                        "(CH3)2CH CH3",     "CH3CH2CH2CH3",    "(CH3)2-CH-CH2-CH3", "CH3(CH2)3CH3", "C6H14" ,
                                        "CH3(CH2)5CH3",     "CCH3-(CH2)6-CH3", "H3C-(CH2)7-CH3",    "CH3(CH2)8CH3", "He" ,      "Ar" };

        std::vector<double> Mm_ = {16.04246, 28.0134, 44.0095, 30.06904, 44.09562, 18.01528, 34.08088, 2.01588, 28.0101, 31.9988, 58.1222, 58.1222, 72.14878, 72.14878, 86.17536, 100.20194, 114.22852, 128.2551, 142.28168, 4.002602, 39.948}; // (g*mol-1)
        std::vector<double> Tc_ = {190.564, 126.192, 304.1282, 305.322, 369.825, 647.096, 373.1, 33.19, 132.86, 154.595, 407.817, 425.125, 460.35, 469.7, 507.82, 540.13, 569.32, 594.55, 617.7, 5.1953, 150.687}; // (K)
        std::vector<double> Pc_ = {4.5992, 3.3958, 7.3773, 4.8722, 4.2512, 22.064, 9, 1.2964, 3.494, 5.043, 3.629, 3.796, 3.378, 3.37, 3.034, 2.736, 2.497, 2.281, 2.103, 0.2276, 4.863}; // (MPa)
        std::vector<double> Acentric_ = {0.01142, 0.0372, 0.22394, 0.0995, 0.1521, 0.3443, 0.1005, -0.219, 0.0497, 0.0222, 0.184, 0.201, 0.2274, 0.251, 0.299, 0.349, 0.393, 0.4433, 0.4884, -0.385, -0.00219};
        std::vector<double> Rhoc_= {10.139342719, 11.1839, 10.624978698, 6.87085454, 5.000043088, 17.87371609, 10.19,    14.94,   10.85,   13.63,   3.86014294, 3.920016792, 3.271,    3.215577588, 2.705877875, 2.315324434, 2.056404127, 1.81,     1.64,      17.399,    13.407429659}; //mol/L

        data dummy;
        dummy.val = 0.0;
        dummy.vec.clear();
        dummy.str = "";

        for (int i =0; i<CAS_.size(); i++){
            dummy.str = CAS_[i];        DBm[std::make_tuple( CAS_[i], GEN, CAS     )]  = dummy;
            dummy.str = name_[i];       DBm[std::make_tuple( CAS_[i], GEN, name    )]  = dummy;
            dummy.str = formula_[i];    DBm[std::make_tuple( CAS_[i], GEN, formula )]  = dummy;
            dummy.str = "";

            dummy.val = Mm_[i];         DBm[std::make_tuple( CAS_[i], GEN, Mm       )] = dummy;
            dummy.val = Tc_[i];         DBm[std::make_tuple( CAS_[i], GEN, Tc       )] = dummy;
            dummy.val = Pc_[i];         DBm[std::make_tuple( CAS_[i], GEN, Pc       )] = dummy;
            dummy.val = Acentric_[i];   DBm[std::make_tuple( CAS_[i], GEN, Acentric )] = dummy;
            dummy.val = Rhoc_[i];       DBm[std::make_tuple( CAS_[i], GEN, Rhoc     )] = dummy;
        }
    }

};


#endif

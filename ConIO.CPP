#include "ConIO.h"
#include "MathUtil.h"

std::map<std::string, std::string> conEsc;

void logHeader(std::ofstream& logs){
    logs << "\n";
    logs << "=====================================================================================\n";
    logs << "========================= N E W   S E S S I O N =====================================\n";
    logs << "=====================================================================================\n";

    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    struct std::tm* ptm = std::localtime(&now);
    logs << "\nStart at: " << std::put_time(ptm, "%c") <<  "\n\n";
    return;
}

void logMix(std::ofstream& logs, std::string name, std::map<std::string, double> mix){
double tot = 0.0;

    logs << "========================== A C T I V E   M I X ======================================\n";
    logs<< "MIXTURE name:          <" << name<<">\n\n";
    for(const auto& gasPair: mix){
        tot += gasPair.second;
        if (gasPair.second!=0.0){
            logs<<std::format("{:20}",gasPair.first)<<" ";
            logs<<std::format("{:7.3f}",gasPair.second)<<'\n';
        }
    }
    logs<<"                       ------\n";
    logs<<"                      "<<std::format("{:7.4f}",tot)<<'\n'<<'\n';
}

void logEOSin(std::ofstream& logs, std::string EOS, double P, double T, double atm){
    logs << "================================= E O S =============================================\n";
    logs << EOS << '\n';
    logs << ">Pressure   = [MPa.g] "<<std::format("{:<9.4f}",P)<<"     [MPa.a] "<<std::format("{:<9.4f}",P+atm)<<"\n";
    logs << ">Pressure   = [psi.g] "<<std::format("{:<9.4f}",P*mpaTOpsi)<<"     [psi.a] "<<std::format("{:<9.4f}",P*mpaTOpsi+atmUS)<<"\n";
    logs << ">Temperature= [C]     "<<std::format("{:<9.3f}",T)<<"     [K]    "<<std::format("{:<9.3f}",Tcon('C',T,'K'))<<"     [F]    "<<std::format("{:<9.3f}",Tcon('C',T,'F'))<<"\n\n";
}

void logEOSdata(std::ofstream& logs, const char& fmt_str, double val, const char& valUM, double Cf=0.0, const char& altUM=' '){
char buffer[1000];

    if(Cf == 0.0) sprintf(buffer, &fmt_str, val, valUM);
    else          sprintf(buffer, &fmt_str, val, valUM, val*Cf, altUM);
}


void setCur(int x = 0 , int y = 0)
{
    HANDLE handle;
    COORD coordinates;
    handle = GetStdHandle(STD_OUTPUT_HANDLE);
    coordinates.X = x;
    coordinates.Y = y;
    SetConsoleCursorPosition ( handle , coordinates );
}

char choice(const char* prompt, const char* valid){
std::string validChar = valid;
std::string text = prompt;
char c;
bool NextHigh;

    fflush(stdin);
    for(char c:text){
        if(c==')') NextHigh= false;
        if (NextHigh) std::cout << "\033[92m"<<c<<"\033[0m";
        else std::cout << c;
        if(c=='(') NextHigh= true;
    }

    do{c=getch();} while (validChar.find_first_of(c)==std::string::npos);
    std::cout << c << '\n';
    return c;
}

int get_int(std::string col)
{
char c;
int ret;
std::string word;

    fflush(stdin);
    std::cout << conEsc[col];
    do {
        c = _getch();
        if (isdigit(c) || (c=='-')|| (c=='+')){
            word.push_back(c);
            std::cout << c;
        }
        if (c=='\b')
            if (!word.empty()){
                word.pop_back();
                std::cout << "\b \b";}
    }while (c != '\r');

    if(word.empty()) {
        ret = 0.0;
        std::cout << ret;
    }
    else ret = stoi(word);
    std::cout << conEsc["Def"];
    return ret;
}

double get_double(std::string col)
{
char c;
double ret;
std::string word;

    fflush(stdin);
    std::cout << conEsc[col];
    do {
        c = _getch();
        if (isdigit(c) || (c=='.')|| (c=='-')|| (c=='+')){
            word.push_back(c);
            std::cout << c;
        }
        if (c=='\b')
            if (!word.empty()){
                word.pop_back();
                std::cout << "\b \b";}
    }while (c != '\r');

    if(word.empty()) {
        ret = 0.0;
        std::cout << ret;
    }
    else ret = stod(word);
    std::cout << conEsc["Def"];
    return ret;
}

std::string get_string(std::string col)
{
char c;
std::string word;

    fflush(stdin);
    std::cout << conEsc[col];
    do {
        c = _getch();
        if (isalpha(c) || isdigit(c)){
            word.push_back(c);
            std::cout << c;
        }
        if (c=='\b')
            if (!word.empty()){
                word.pop_back();
                std::cout << "\b \b";}
    }while (c != '\r');
    std::cout << conEsc["Def"];
    return word;
}

char conMenu(const char* Title, std::vector<std::string>& menuList){
int scrx, scry;
bool NextHigh=false;
std::string validChar;
char c;

    scrx = 4;
    scry = 3;
    system("cls");
    std::cout <<"================================================================================\n";
    std::cout <<Title;
    for(std::string& item:menuList){
        setCur(scrx,scry++);
        for(char c:item){
            if (NextHigh){
                std::cout << "\033[92m"<<c<<"\033[0m";
                validChar.push_back(c);
                NextHigh= false;
            }
            else std::cout << c;
            if(c=='(') NextHigh= true;
        }
        std::cout << '\n';
    }
    fflush(stdin);
    setCur(scrx,++scry);
    std::cout << "Select -> ";
    do{c=getch();} while (validChar.find_first_of(c)==std::string::npos);
    std::cout << "\033[92m" << c << "\033[0m" <<'\n';
    return c;
}

void initCon(){
conEsc["BlackF"] =       "\033[30m";
conEsc["BlackB"] =       "\033[40m";
conEsc["RedF"] =         "\033[31m";
conEsc["RedB"] =         "\033[41m";
conEsc["GreenF"] =       "\033[32m";
conEsc["GreenB"] =       "\033[42m";
conEsc["YellowF"] =      "\033[33m";
conEsc["YellowB"] =      "\033[43m";
conEsc["BlueF"] =        "\033[34m";
conEsc["BlueB"] =        "\033[44m";
conEsc["MagentaF"] =     "\033[35m";
conEsc["MagentaB"] =     "\033[45m";
conEsc["CyanF"] =        "\033[36m";
conEsc["CyanB"] =        "\033[46m";
conEsc["LhGreenF"] =     "\033[37m"; //LightGreenF
conEsc["LhGreenB"] =     "\033[47m";
conEsc["DrGrayF"] =      "\033[90m"; //Dark gray
conEsc["DrGrayB"] =      "\033[100m";
conEsc["LhRedF"] =       "\033[91m"; //LightRed
conEsc["LhRedB"] =       "\033[101m";
conEsc["LhGreenF"] =     "\033[92m";
conEsc["LhGreenB"] =     "\033[102m";
conEsc["LhYellowF"] =    "\033[93m";
conEsc["LhYellowB"] =    "\033[103m";
conEsc["LhBlueF"] =      "\033[94m";
conEsc["LhBlueB"] =      "\033[104m";
conEsc["LhMagentaF"] =   "\033[95m";
conEsc["LhMagentaB"] =   "\033[105m";
conEsc["LhCyanF"] =      "\033[96m";
conEsc["LhCyanB"] =      "\033[106m";
conEsc["WhiteF"] =       "\033[97m";
conEsc["WhiteB"] =       "\033[107m";
conEsc["Bold"] =         "\033[1m";
conEsc["Underline"] =    "\033[4m";
conEsc["NoUnderline"] =  "\033[24m";
conEsc["Negative"] =     "\033[7m";  //(reverse the foreground and background)
conEsc["Positive"] =     "\033[27m"; //(no negative)
conEsc["Def"] =          "\033[0m"; //Default
}

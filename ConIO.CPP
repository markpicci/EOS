#include "ConIO.h"
#include "MathUtil.h"

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

void verbose(std::ofstream& logs, std::string print, bool cr){
    if(cr){
        logs << print << "\n";
        std::cout << print << "\n";
    } else {
        logs << print;
        std::cout << print;
    }
}

void setCurPos(int x, int y)
{
    COORD coordinates = {(SHORT)x, (SHORT)y};
    HANDLE handle;
    handle = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleCursorPosition ( handle , coordinates );
}

void setCurPos(COORD coordinates)
{
    HANDLE handle;
    handle = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleCursorPosition (handle, coordinates);
}

COORD getCurPos()
{
    HANDLE handle;
    handle = GetStdHandle(STD_OUTPUT_HANDLE);

    CONSOLE_SCREEN_BUFFER_INFO cbsi;
    if (GetConsoleScreenBufferInfo(handle, &cbsi))
    {
        return cbsi.dwCursorPosition;
    }
    else
    {
        // The function failed. Call GetLastError() for details.
        COORD invalid = {0, 0 };
        return invalid;
    }
}

void clearScreen(char fill) {
    HANDLE handle;
    handle = GetStdHandle(STD_OUTPUT_HANDLE);

    COORD Orig = {0,0};
    CONSOLE_SCREEN_BUFFER_INFO s;
    GetConsoleScreenBufferInfo(handle, &s);
    DWORD written, cells = s.dwSize.X * s.dwSize.Y;
    FillConsoleOutputCharacter(handle, fill, cells, Orig, &written);
    FillConsoleOutputAttribute(handle, s.wAttributes, cells, Orig, &written);
    SetConsoleCursorPosition(handle, Orig);
}

void MaximizeWindow()
{
	CONSOLE_SCREEN_BUFFER_INFO info;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);
	SMALL_RECT rc;
	rc.Left = rc.Top = 0;
	rc.Right = (short)(std::min(info.dwMaximumWindowSize.X, info.dwSize.X) - 1);
	rc.Bottom = (short)(std::min(info.dwMaximumWindowSize.Y, info.dwSize.Y) - 1);
	SetConsoleWindowInfo(GetStdHandle(STD_OUTPUT_HANDLE), true, &rc);
}

char choice(const std::string text, std::string validChar, char defChar){
char c;
bool NextHigh = false, first = true;

    COORD coordinates = getCurPos();
    for(size_t i=0; i<text.size()+1; i++) std::cout << " ";
    setCurPos(coordinates);

    fflush(stdin);
    for(char c:text){
        if(c==')') NextHigh= false;
        if (NextHigh) std::cout << "\033[92m"<<c<<"\033[0m";
        else std::cout << c;
        if(c=='(') NextHigh= true;
    }

    std::cout << "\033[92m";
    if (defChar!=0) std::cout << defChar << "\b";
    do{
        c=getch();
        if (first && defChar!=0 && c=='\r')
            {std::cout << "\033[0m\n"; return defChar;}
        first=false;
    } while (validChar.find_first_of(c)==std::string::npos);
    std::cout << c << "\033[0m\n";
    return c;
}

std::string getStrGen(std::string validChar, std::string col,std::string defStr = "", size_t charLen=20)
{
char c;
int ret;
std::string word;

    //std::cout << conEsc["DrGrayB"];
    std::cout << getEsc(col);
    COORD coordinates = getCurPos();
    for(size_t i=0; i<charLen; i++) std::cout << "_";
    setCurPos(coordinates);
    if(!defStr.empty()) {
        word = defStr.substr(0, charLen);
        std::cout << word;
    }

    fflush(stdin);
    do {
        c = _getch();
        if (validChar.find_first_of(c)!=std::string::npos && word.size()<charLen){
            word.push_back(c);
            std::cout << c;
        }
        if (c=='\b')
            if (!word.empty()){
                word.pop_back();
                std::cout << "\b_\b";}
    }while (c != '\r');

    return word;
}

int get_int(std::string col, int defInt)
{
char c;
int ret;
std::string word;

    if (defInt == 0) word = ""; else word = std::format("{:d}", defInt);
    word=getStrGen("01234567890-+", col, word);

    if(word.empty()) {
        ret = 0.0;
        std::cout << ret;
    }
    else ret = stoi(word);
    std::cout << getEsc("Def");
    return ret;
}

double get_double(std::string col, double defDouble)
{
char c;
double ret, fracPart, intPart;
std::string word;


    if (defDouble==0.0){
        word = "";
    } else {
        if (abs(defDouble)<0.00001)
             word = std::format("{:.16f}", defDouble);
        else word = std::format("{:f}", defDouble);
        word.erase(std::find_if_not(word.rbegin(), word.rend(), [](char c){ return c == '0'; }).base(), word.end());
        if (word.at(word.size()-1)=='.') word+="0";
    }

    word=getStrGen("01234567890.-+", col, word);

    if(word.empty()) {
        ret = 0.0;
        std::cout << ret;
    }
    else ret = stod(word);
    std::cout << getEsc("Def");
    return ret;
}

std::string get_string(std::string col, std::string defStr)
{
char c;
std::string word;

    word=getStrGen("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz01234567890.-+_", col, defStr);

    std::cout << getEsc("Def");
    return word;
}

char conMenu(const char* Title, std::vector<std::string>& menuList){
int scrx, scry;
bool NextHigh=false;
std::string validChar;
char c;

    scrx = 4;
    scry = 3;
    clearScreen();
    std::cout <<"================================================================================\n";
    std::cout <<Title;
    for(std::string& item:menuList){
        setCurPos(scrx,scry++);
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
    setCurPos(scrx,++scry);
    std::cout << "Select -> ";
    do{c=getch();} while (validChar.find_first_of(c)==std::string::npos);
    std::cout << "\033[92m" << c << "\033[0m" <<'\n';
    return c;
}

std::string getEsc(const std::string& col){
    if(col=="BlackF")       return "\033[30m";
    if(col=="BlackB")       return "\033[40m";
    if(col=="RedF")         return "\033[31m";
    if(col=="RedB")         return "\033[41m";
    if(col=="GreenF")       return "\033[32m";
    if(col=="GreenB")       return "\033[42m";
    if(col=="YellowF")      return "\033[33m";
    if(col=="YellowB")      return "\033[43m";
    if(col=="BlueF")        return "\033[34m";
    if(col=="BlueB")        return "\033[44m";
    if(col=="MagentaF")     return "\033[35m";
    if(col=="MagentaB")     return "\033[45m";
    if(col=="CyanF")        return "\033[36m";
    if(col=="CyanB")        return "\033[46m";
    if(col=="LhGreenF")     return "\033[37m";  //LightGreenF
    if(col=="LhGreenB")     return "\033[47m";
    if(col=="DrGrayF")      return "\033[90m";  //Dark gray
    if(col=="DrGrayB")      return "\033[100m";
    if(col=="LhGrayF")      return "\033[37m";  //Ligt gray
    if(col=="LhGrayB")      return "\033[47m";
    if(col=="LhRedF")       return "\033[91m";  //LightRed
    if(col=="LhRedB")       return "\033[101m";
    if(col=="LhGreenF")     return "\033[92m";
    if(col=="LhGreenB")     return "\033[102m";
    if(col=="LhYellowF")    return "\033[93m";
    if(col=="LhYellowB")    return "\033[103m";
    if(col=="LhBlueF")      return "\033[94m";
    if(col=="LhBlueB")      return "\033[104m";
    if(col=="LhMagentaF")   return "\033[95m";
    if(col=="LhMagentaB")   return "\033[105m";
    if(col=="LhCyanF")      return "\033[96m";
    if(col=="LhCyanB")      return "\033[106m";
    if(col=="WhiteF")       return "\033[97m";
    if(col=="WhiteB")       return "\033[107m";
    if(col=="Bold")         return "\033[1m";
    if(col=="Underline")    return "\033[4m";
    if(col=="NoUnderline")  return "\033[24m";
    if(col=="Negative")     return "\033[7m";   //(reverse the foreground and background)
    if(col=="Positive")     return "\033[27m";  //(no negative)
    return "\033[0m";                           //Default
}

// ============================================================================
// Color methods
// ============================================================================

WORD GetTextColor(Color color){
    //CONSOLE_SCREEN_BUFFER_INFO csbi;
    //GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    //WORD attributes = csbi.wAttributes & ~FOREGROUND_RED & ~FOREGROUND_GREEN & ~FOREGROUND_BLUE & ~FOREGROUND_INTENSITY;
    WORD attributes = 0;

    switch (color)
    {
        case BLACK:           attributes |= 0;                                                                            break;
        case LIGHT_RED:       attributes |= FOREGROUND_RED   | FOREGROUND_INTENSITY;                                      break;
        case LIGHT_GREEN:     attributes |= FOREGROUND_GREEN | FOREGROUND_INTENSITY;                                      break;
        case LIGHT_BLUE:      attributes |= FOREGROUND_BLUE  | FOREGROUND_INTENSITY;                                      break;
        case LIGHT_CYAN:      attributes |= FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY;                    break;
        case LIGHT_MAGENTA:   attributes |= FOREGROUND_RED   | FOREGROUND_BLUE | FOREGROUND_INTENSITY;                    break;
        case LIGHT_YELLOW:    attributes |= FOREGROUND_GREEN | FOREGROUND_RED  | FOREGROUND_INTENSITY;                    break;
        case LIGHT_GRAY:      attributes |= FOREGROUND_RED   | FOREGROUND_GREEN | FOREGROUND_BLUE;                        break;
        case DARK_RED:        attributes |= FOREGROUND_RED;                                                               break;
        case DARK_GREEN:      attributes |= FOREGROUND_GREEN;                                                             break;
        case DARK_BLUE:       attributes |= FOREGROUND_BLUE;                                                              break;
        case DARK_CYAN:       attributes |= FOREGROUND_GREEN | FOREGROUND_BLUE;                                           break;
        case DARK_MAGENTA:    attributes |= FOREGROUND_RED   | FOREGROUND_BLUE;                                           break;
        case DARK_YELLOW:     attributes |= FOREGROUND_GREEN | FOREGROUND_RED;                                            break;
        case DARK_GRAY:       attributes |= FOREGROUND_INTENSITY;                                                         break;
        case WHITE:           attributes |= FOREGROUND_RED   | FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY; break;
    }

    return attributes;
}

void SetTextColor(Color color) {
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), GetTextColor(color));
}

WORD GetBackgroundColor(Color color) {
    //CONSOLE_SCREEN_BUFFER_INFO csbi;
    //GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE);, &csbi);
    //WORD attributes = csbi.wAttributes & ~BACKGROUND_RED & ~BACKGROUND_GREEN & ~BACKGROUND_BLUE & ~BACKGROUND_INTENSITY;
    WORD attributes = 0;

    switch (color)
    {
        case BLACK:           attributes |= 0;                                                                            break;
        case LIGHT_RED:       attributes |= BACKGROUND_RED   | BACKGROUND_INTENSITY;                                      break;
        case LIGHT_GREEN:     attributes |= BACKGROUND_GREEN | BACKGROUND_INTENSITY;                                      break;
        case LIGHT_BLUE:      attributes |= BACKGROUND_BLUE  | BACKGROUND_INTENSITY;                                      break;
        case LIGHT_CYAN:      attributes |= BACKGROUND_GREEN | BACKGROUND_BLUE  | BACKGROUND_INTENSITY;                   break;
        case LIGHT_MAGENTA:   attributes |= BACKGROUND_RED   | BACKGROUND_BLUE  | BACKGROUND_INTENSITY;                   break;
        case LIGHT_YELLOW:    attributes |= BACKGROUND_GREEN | BACKGROUND_RED   | BACKGROUND_INTENSITY;                   break;
        case LIGHT_GRAY:      attributes |= BACKGROUND_RED   | BACKGROUND_GREEN | BACKGROUND_BLUE;                        break;
        case DARK_RED:        attributes |= BACKGROUND_RED;                                                               break;
        case DARK_GREEN:      attributes |= BACKGROUND_GREEN;                                                             break;
        case DARK_BLUE:       attributes |= BACKGROUND_BLUE;                                                              break;
        case DARK_CYAN:       attributes |= BACKGROUND_GREEN | BACKGROUND_BLUE;                                           break;
        case DARK_MAGENTA:    attributes |= BACKGROUND_RED   | BACKGROUND_BLUE;                                           break;
        case DARK_YELLOW:     attributes |= BACKGROUND_GREEN | BACKGROUND_RED;                                            break;
        case DARK_GRAY:       attributes |= BACKGROUND_INTENSITY;                                                         break;
        case WHITE:           attributes |= BACKGROUND_RED   | BACKGROUND_GREEN | BACKGROUND_BLUE | BACKGROUND_INTENSITY; break;
    }

    return attributes;
}

void SetBackgroundColor(Color color) {
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), GetBackgroundColor(color));
}

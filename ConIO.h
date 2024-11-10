#ifndef __CON_IO__H
#define __CON_IO__H

#include <map>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <conio.h>
#include <string>
#include <vector>
#include <chrono>
#include <windows.h>
#include <string_view>
#include <format>

enum Color{
    BLACK,
    LIGHT_RED,
    LIGHT_GREEN,
    LIGHT_BLUE,
    LIGHT_CYAN,
    LIGHT_MAGENTA,
    LIGHT_YELLOW,
    LIGHT_GRAY,
    DARK_RED,
    DARK_GREEN,
    DARK_BLUE,
    DARK_CYAN,
    DARK_MAGENTA,
    DARK_YELLOW,
    DARK_GRAY,
    WHITE
};

char        conMenu(const char* Title, std::vector<std::string>& menuList);
std::string getEsc(const std::string& col);
char        choice(const std::string text, std::string validChar, char defChar=0);
int         get_int   (std::string col = "WhiteF", int defInt = 0);
double      get_double(std::string col = "WhiteF", double defDouble = 0.0);
std::string get_string(std::string col = "WhiteF", std::string defStr = "");

void        logHeader(std::ofstream& logs);
void        logMix(std::ofstream& logs, std::string name, std::map<std::string, double> mix);
void        logEOSin(std::ofstream& logs, std::string EOS, double P, double T, double atm);
void        verbose(std::ofstream& logs, std::string print, bool cr = false);

COORD getCurPos();
void  setCurPos(int x = 0 , int y = 0);
void  setCurPos(COORD coordinates);
WORD  GetTextColor(Color color);
void  SetTextColor(Color color);
WORD  GetBackgroundColor(Color color);
void  SetBackgroundColor(Color color);
void  clearScreen(char fill = ' ');
void  MaximizeWindow();

#endif


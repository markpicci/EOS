#ifndef __CON-IO__
#define __CON-IO__

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

void initCon(void);
int get_int(std::string col = "WhiteF");
double get_double(std::string col = "WhiteF");
std::string get_string(std::string col = "WhiteF");
char conMenu(const char* Title, std::vector<std::string>& menuList);
char choice(const char* prompt, const char* valid);
void logHeader(std::ofstream& logs);
void logMix(std::ofstream& logs, std::string name, std::map<std::string, double> mix);
void logEOSin(std::ofstream& logs, std::string EOS, double P, double T, double atm);
#endif


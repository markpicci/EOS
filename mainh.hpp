#ifndef __MAINH__
#define __MAINH__

#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <stdlib.h>
#include <conio.h>
#include <filesystem>
#include <vector>

struct mixSt{
    std::string name;
    std::map<std::string, double> mix;
};

void equipmentMain(std::ofstream& logs);
void EOS(std::ofstream& logs, const size_t activeMix, std::vector<mixSt>& mix);

#endif

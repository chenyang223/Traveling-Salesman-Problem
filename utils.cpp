#include "utils.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;
std::string &GetDataFile()
{
    static std::string data_file;
    return data_file;
}

std::string &GetDataPath()
{
    static std::string data_path;
    return data_path;
}

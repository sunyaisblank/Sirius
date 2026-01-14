// CRCF001A.h - Configuration Class

#pragma once

#include <string>
#include <vector>
#include <map>

struct Param {
    double value;
    double min;
    double max;
};

using Config = std::map<std::string, Param>;
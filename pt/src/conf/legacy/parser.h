#pragma once
#include "legacy.pb.h"
#include <istream>

namespace conf::legacy {
    class Parser {
    public:
        LegacyConf parse(std::istream& is);
    };
}
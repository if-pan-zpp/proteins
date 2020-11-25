#pragma once
#include "legacy.pb.h"
#include <istream>

namespace conf {
    class LegacyParser {
    public:
        LegacyConf parse(std::istream& is);
    };
}
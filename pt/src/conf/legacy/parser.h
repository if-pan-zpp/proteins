#pragma once
#include "legacy.pb.h"
#include <istream>

namespace conf::legacy {
    class Parser {
    public:
        Conf parse(std::istream& is);
    };
}
#pragma once
#include "ir.pb.h"
#include <string>

namespace pdb::ir {
    class LineParser {
    public:
        virtual void parse(std::string_view line, Field *field) = 0;
    };
}

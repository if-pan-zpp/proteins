#pragma once
#include "ir.pb.h"
#include <string>

namespace pdb::ir {
    class LineParser {
    public:
        virtual void parse(std::string_view line, AnyField *field) = 0;
    };
}

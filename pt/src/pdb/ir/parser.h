#pragma once
#include <istream>
#include <string_view>
#include <memory>
#include <unordered_map>
#include "ir.pb.h"
#include "line_parser.h"

namespace pdb::ir {
    class Parser {
    private:
        std::unordered_map<Headers, std::unique_ptr<LineParser>> subparsers;

    public:
        Parser();
        Entry parse(std::istream& is);
    };
}

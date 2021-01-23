#include "parser.hpp"
#include "auto_line_parser.hpp"
#include "util/str.hpp"
#include <string>
#include <string_view>
using namespace std;
using namespace pdb::ir;

Parser::Parser() {
    subparsers[Headers::ATOM] = make_unique<AutoLineParser>(
        Headers::ATOM,
        AutoLineParser::Sectors {
            {7, 11}, {13, 16}, {17, 17}, {18, 20}, {22, 22}, {23, 26},
            {27, 27}, {31, 38}, {39, 46}, {47, 54}, {55, 60}, {61, 66},
            {77, 78}, {79, 80}
        });

    subparsers[Headers::TER] = make_unique<AutoLineParser>(
        Headers::TER,
        AutoLineParser::Sectors {
            {7, 11}, {18, 20}, {22, 22}, {23, 26}, {27, 27}
        });

    subparsers[Headers::CRYST1] = make_unique<AutoLineParser>(
        Headers::CRYST1,
        AutoLineParser::Sectors {
            {7, 15}, {16, 24}, {25, 33}, {34, 40}, {41, 47}, {48, 54},
            {56, 66}, {67, 70}
        });

    subparsers[Headers::SSBOND] = make_unique<AutoLineParser>(
        Headers::SSBOND,
        AutoLineParser::Sectors {
            {8, 10}, {12, 14}, {16, 16}, {18, 21}, {22, 22}, {26, 28},
            {30, 30}, {32, 35}, {46, 46}, {60, 65}, {67, 72}, {74, 78}
        });

    subparsers[Headers::END] = make_unique<AutoLineParser>(
        Headers::END,
        AutoLineParser::Sectors {
        });
}

Entry Parser::parse(std::istream &is) {
    Entry entry;
    Field* field;
    string line;

    while (getline(is, line)) {
        string header_name(line, 0, 6);
        header_name = util::rtrim(header_name);

        auto header_type = Headers_descriptor()->FindValueByName(header_name);
        if (!header_type) continue;
        auto header = (Headers)header_type->number();

        field = entry.add_fields();
        subparsers.at(header)->parse(line, field);
        if (field->has_end()) break;
    }

    return entry;
}

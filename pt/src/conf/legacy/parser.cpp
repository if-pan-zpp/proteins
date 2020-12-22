#include "parser.h"
#include "util/str.h"
#include "util/field.h"
#include <string>
using namespace conf::legacy;
using namespace std;
using namespace google::protobuf;

Conf Parser::parse(istream &is) {
    string line;
    Conf conf;
    auto refl = Conf::GetReflection();
    auto desc = Conf::GetDescriptor();

    while (getline(is, line)) {
        auto parts = ::util::split(line);
        if (parts.size() < 2) continue;

        auto name = string(parts[0]);
        auto field = desc->FindFieldByName(name);
        if (!field) continue;

        auto ufield = ::util::Field(&conf, field, refl);
        ufield.setFrom(parts[1]);
    }

    return conf;
}

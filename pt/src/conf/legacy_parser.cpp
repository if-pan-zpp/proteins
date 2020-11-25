#include "legacy_parser.h"
#include "util/strops.h"
#include "util/field.h"
#include <string>
using namespace conf;
using namespace std;
using namespace google::protobuf;

LegacyConf LegacyParser::parse(istream &is) {
    string line;
    LegacyConf conf;
    auto refl = conf.GetReflection();

    while (getline(is, line)) {
        auto parts = ::util::split(line);
        if (parts.size() < 2) continue;

        auto name = string(parts[0]);
        auto field = conf.descriptor()->FindFieldByName(name);
        if (!field) continue;

        auto ufield = ::util::Field(&conf, field, refl);
        ufield.setFrom(parts[1]);
    }

    return conf;
}

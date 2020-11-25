#include <fstream>
#include "conf/legacy_parser.h"
using namespace std;
using namespace pdb::ir;

int main() {
    ifstream config("inputfile1");
    conf::LegacyParser parser;
    parser.parse(config).PrintDebugString();
    return 0;
}

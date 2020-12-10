#include <fstream>
#include "pdb/ir/parser.h"
using namespace std;

int main() {
    ifstream config("1ubq.pdb");
    pdb::ir::Parser parser;
    auto ir = parser.parse(config);
    ir.PrintDebugString();
    return 0;
}

#include <fstream>
#include "pdb/ir/ir_parser.h"
using namespace std;
using namespace pdb::ir;

int main() {
    ifstream ubq("1ubq.pdb");
    auto irp = IRParser();
    auto ir = irp.parse(ubq);
    return 0;
}

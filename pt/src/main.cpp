#include <fstream>
#include "conf/legacy/parser.h"
using namespace std;

int main() {
    ifstream config("inputfile1");
    conf::legacy::Parser parser;
    parser.parse(config).PrintDebugString();
    return 0;
}

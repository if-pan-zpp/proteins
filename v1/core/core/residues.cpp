#include "residues.h"
#include <unordered_map>
#include <algorithm>
using namespace std;

Residue::Residue(char code) {
    this->code = code;
}

Residue::Residue(std::string name) {
    static std::unordered_map<std::string, char> nameToCode {
        {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'}, {"CYS", 'C'},
        {"GLU", 'E'}, {"GLN", 'Q'}, {"GLY", 'G'}, {"HIS", 'H'}, {"ILE", 'I'},
        {"LEU", 'L'}, {"LYS", 'K'}, {"MET", 'M'}, {"PHE", 'F'}, {"PRO", 'P'},
        {"SER", 'S'}, {"THR", 'T'}, {"TRP", 'W'}, {"TYR", 'Y'}, {"VAL", 'V'}
    };

    transform(name.begin(), name.end(), name.begin(), ::toupper);
    code = nameToCode[name];
}

inline Residue::operator char() const {
    return code;
}

Residue::operator std::string const&() const {
    static std::unordered_map<char, std::string> codeToName {
        {'A', "ALA"}, {'R', "ARG"}, {'N', "ASN"}, {'D', "ASP"}, {'C', "CYS"},
        {'E', "GLU"}, {'Q', "GLN"}, {'G', "GLY"}, {'H', "HIS"}, {'I', "ILE"},
        {'L', "LEU"}, {'K', "LYS"}, {'M', "MET"}, {'F', "PHE"}, {'P', "PRO"},
        {'S', "SER"}, {'T', "THR"}, {'W', "TRP"}, {'Y', "TYR"}, {'V', "VAL"}
    };

    return codeToName[code];
}

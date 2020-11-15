#pragma once
#include <string>
#include <unordered_map>

class Residue {
private:
    char code;

public:
    Residue(char code);
    Residue(std::string name);

    inline operator char() const;
    operator std::string const&() const;
};

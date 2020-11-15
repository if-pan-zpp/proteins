#pragma once
#include <istream>
#include <ostream>
#include <string>
#include <Eigen/Core>

class Atom {
public:
    Eigen::Vector3d r;
    std::string name;
    double occupany, tempFactor;
};

class Residue {
    
};

class PDB {
public:
    friend std::istream& operator>>(std::istream& is, PDB& pdb);
    friend std::ostream& operator<<(std::ostream& os, PDB const& pdb);
};

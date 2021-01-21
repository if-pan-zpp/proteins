#include <fstream>
#include "pdb/ir/parser.h"
#include "state/def.h"
#include "conf/constants.h"
#include <eigen3/Eigen/Core>
using namespace std;

class Config {
public:
    Config(int argc, char **argv);

    // Angle potential:
    bool angle_pot = true;
    bool dih_angle_pot = false;
    bool tabularized_angle_pot;

    // Other configs:
    string out_filename;
    // ...
};

class State {
public:
    
    void init(const Config &config);
};

class VerList {
public:
    
    void init(const Config &config);
};

/* Class for storing pseudo-atoms.
 */
class PAtoms {
public:
    void init(const Config &config);
    int n;
    array<Vec3DArray, DER_ORDER + 1> der; // positions' derivatives.
    
};

class Statistics {
public:
    void init(const Config &config);
};

class Integrator {
public:
    void init(const Config &config, const PAtoms &p_atoms);
    void step(PAtoms &p_atoms);
private:
};

// POTENTIALS

class Potential {
public:
    virtual Vec3DArray calculate_forces(
	const State &state,
	const PAtoms &p_atoms,
	const VerList &verlet_list,
	Statistics &stats) = 0;

    virtual bool is_enabled(const Config &config) const = 0;
    virtual void init(const Config &config) {};
};

class BondAnglePotential : public Potential {
public:
    Vec3DArray calculate_forces(
	const State &state,
	const PAtoms &p_atoms,
	const VerList &verlet_list,
	Statistics &stats) override { return Vec3DArray(); }
    bool is_enabled(const Config &config) const override; 
    void init(const Config &config) override {}
};

bool BondAnglePotential::is_enabled(const Config &config) const {
    return config.angle_pot;
}

class DihAnglePotential : public Potential {
public:
    Vec3DArray calculate_forces(
	const State &state,
	const PAtoms &p_atoms,
	const VerList &verlet_list,
	Statistics &stats) override { return Vec3DArray(); }
    bool is_enabled(const Config &config) const override; 
    void init(const Config &config) override {}
};

bool DihAnglePotential::is_enabled(const Config &config) const {
    return config.dih_angle_pot;
}

vector<unique_ptr<Potential>> all_supported_potentials() {
    vector<unique_ptr<Potential>> all;
    all.emplace_back(make_unique<BondAnglePotential>());
    all.emplace_back(make_unique<DihAnglePotential>());
    return all;
}


/*
  This class deals with everything that
  isn't taking/checking some input from a user.
 */
class Simulation {
public:
    Simulation(const Config &config);
    void run();
private:

    State state;
    VerList verlet_list;    
    PAtoms p_atoms; 
    vector<shared_ptr<Potential>> potentials;
    Integrator integrator;
};

Simulation::Simulation(const Config &config) {
    for (auto &pot_ptr: all_supported_potentials()) {
	if (pot_ptr -> is_enabled(config)) potentials.emplace_back(move(pot_ptr));
    }

    for (auto &pot_ptr : potentials) pot_ptr -> init(config);

    cout << "Number of potentials: " << potentials.size() << endl;
}

void Simulation::run() {}

Config::Config(int argc, char **argv) {}

int main(int argc, char **argv) {
    Config config(argc, argv);
    Simulation simulation(config);
    simulation.run();
    return 0;
}

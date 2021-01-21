/*
  Janek:
  Integrator class - easy

  Derived classes of potentials
    - harmonic between (i,i+1)
    - repulsive (i,i+2)

  Kuba Bo:
    - bond angle
    - dihedral angle
    - repulsive any pair <- uses verlet list
    - L-J between native contacts 

  Paweł:
  Reading PDB and configs and getting native conformation, amino acid sequence etc.

  Kuba Be:
  Verlet List - easy

  


 */


#include <fstream>
#include "pdb/ir/parser.h"
#include "state/def.h"
#include "def/types.h"
#include "conf/constants.h"
#include <eigen3/Eigen/Core>
using namespace std;

/*
 * This should cover all configuration that isn't done at compile time. 
 * It should always be constant after construction.
 *
 * Almost every class will get a const ref to this class in its constructor.
 * But it is discouraged to keep this ref, rather we should copy the relevant
 * parameters, because in most cases it will be a handful of variables and this
 * way we clearly show in one place which parameters influence this particular class.
 */
class Config {
public:
    Config(int argc, char **argv) {
	// Here there will be a lot of parsing configs, checking them,
	// printing to the user, prompting the user...
    }

    // Angle potential:
    bool angle_pot = true;
    bool dih_angle_pot = false;
    bool tabularized_angle_pot;

    // Other configs:
    string out_filename;
    int max_steps = 10;
    def::Scalar verlet_list_max_eps = 10.0;
    // ...
};

/*
 * Class that stores the state of the environment. 
 * Examples: temperature, movement of walls, current step.
 *
 * We need to figure out which classes should be allowed to make changes here.
 */
class State {
public:
    State(const Config &config) {
	max_steps = config.max_steps;
    }
    bool is_running() {
	return cur_step < max_steps;
    }
    void take_step() {
	cur_step++;
    }
    int max_steps = 0;
    int cur_step = 0;
    def::Scalar temperature = 0.15;
};


/*
 * Class for storing pseudo-atoms. It holds:
 * - positions and their derivatives
 * - masses, amino acid sequence etc.
 * - division into chains
 * - division into structured and unstructured parts
 * - method for calculating the distance between two pseudatoms
 *   (because it depends on periodic boundary conditions and is used a lot)
 * - ...
 * 
 * After construction, only Integrator can make changes here.
 */
class PAtoms {
public:
    PAtoms(const Config &config) {
	n = 1;
    }
    int n;
    array<Vec3DArray, DER_ORDER + 1> der; // positions' derivatives.
    
};

/*
 * Class that provides list of pairs of pseudoatoms that are close.
 * The exact design and implementation is a very hard and interesting problem.
 * But for the purposes of the prototype that follows cg.f in its implementation,
 * it's easy for now.
 */
class VerList {
public:
    VerList(const Config &_config, const PAtoms &_p_atoms)
	: p_atoms(_p_atoms) {

	eps_upper_bound = _config.verlet_list_max_eps;
    }

    using VerIt = vector<pair<int,int>>::const_iterator;

    /* 
     * This returns a range in an unmodifiable vector of pairs
     * that contains every pair closer than 'eps'. 
     */
    pair<VerIt, VerIt> get_verlet_list(def::Scalar eps) const {
	return {list.begin(), list.end()};
    }

    /* 
     * This calculates the distances between last saved positions and current positions.
     * If it's too big, it recalculates the list.
    */
    void take_step() {
    }

    /* 
     * Potential objects can register the eps that they need before the start of the simulation.
     */
    void register_eps(def::Scalar req_eps) {
	assert (req_eps < eps_upper_bound);
	biggest_req_eps = max(biggest_req_eps, req_eps);
    }

private:
    const PAtoms &p_atoms;
    def::Scalar eps_upper_bound;
    def::Scalar biggest_req_eps = 0.;
    vector<pair<int, int>> list;
  
};


/* 
 * 
 */
class Statistics {
public:
    Statistics(const Config &_config,
	       const State &_state,
	       const PAtoms &_p_atoms)
	: state(_state), p_atoms(_p_atoms) {};

    void take_step() {
	cout << "We're in step " << state.cur_step << endl;
	
    }

private:
    const State &state;
    const PAtoms &p_atoms;
};

/*
 * This deals with langevin dynamics (damping and thermal noise)
 * as well as the predictor-corrector algorithm.
 */
class Integrator {
public:
    Integrator(const Config &_config,
	       const State &_state,
	       const PAtoms &_p_atoms) : state(_state) {}

    void take_step(PAtoms &p_atoms, const Vec3DArray &forces) {
	def::Scalar temp = state.temperature;
    }

private:
    const State &state;
};

// POTENTIALS


class Potential {
public:
    Potential(const Config &_config,
	      const State &_state,
	      const PAtoms &_p_atoms) :
	state(_state), p_atoms(_p_atoms) {}

    virtual void init_and_register(VerList &verlet_list) {}

    /*
     * We split one step into those three functions, because
     * in the future we may want to be able to write calculate_forces
     * in a way that can run concurrently with other calculate_forces calls
     * (from this and other potentials). But some elements need synchronization anyway. 
     */
    virtual void init_step() {}
    virtual Vec3DArray calculate_forces(const VerList &verlet_list) = 0;
    virtual void finish_step(Statistics &statistics) {}

    virtual bool is_enabled() const = 0;

protected:
    const State &state;
    const PAtoms &p_atoms;
};

class BondAnglePotential : public Potential {
public:
    BondAnglePotential(const Config &_config,
		       const State &_state,
		       const PAtoms &_p_atoms) :
	Potential(_config, _state, _p_atoms) {

	enabled = _config.angle_pot;
    }

    Vec3DArray calculate_forces(const VerList &verlet_list) override {
	return Vec3DArray::Zero(p_atoms.n, 3);
    }
    bool is_enabled() const override {
	return enabled;
    }
private:
    bool enabled = false;
};

class DihAnglePotential : public Potential {
public:
    DihAnglePotential(const Config &_config,
		      const State &_state,
		      const PAtoms &_p_atoms) :
	Potential(_config, _state, _p_atoms) {

	enabled = _config.dih_angle_pot;
    }

    void init_and_register(VerList &verlet_list) override {
	// Right now, if this potential is enabled, this will fail,
	// because maximum eps is 10.0.
	verlet_list.register_eps(15.0);
    }
    Vec3DArray calculate_forces(const VerList &verlet_list) override {
	return Vec3DArray::Zero(p_atoms.n, 3);
    }
    bool is_enabled() const override {
	return enabled;
    }
private:
    bool enabled = false;
};

// END OF POTENTIALS



/*
 * This class deals with everything that
 * isn't taking/checking some input from a user.
 * In particular, it runs the main loop.
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
    Statistics statistics;

    Vec3DArray forces;

    void take_step();
    vector<unique_ptr<Potential>> all_supported_potentials(const Config &config) const;
};

vector<unique_ptr<Potential>> Simulation::all_supported_potentials(const Config &config) const {
    vector<unique_ptr<Potential>> all;
    all.emplace_back(make_unique<BondAnglePotential>(config, state, p_atoms));
    all.emplace_back(make_unique<DihAnglePotential>(config, state, p_atoms));
    return all;
}

Simulation::Simulation(const Config &config) :
    state(config),
    p_atoms(config),
    verlet_list(config, p_atoms),
    integrator(config, state, p_atoms),
    statistics(config, state, p_atoms)
    {
    
    for (auto &pot_ptr: all_supported_potentials(config)) {
	if (pot_ptr -> is_enabled()) potentials.emplace_back(move(pot_ptr));
    }

    for (auto &pot_ptr : potentials) pot_ptr -> init_and_register(verlet_list);

    cout << "Number of enabled potentials: " << potentials.size() << endl;
}

void Simulation::take_step() {
    state.take_step();
    verlet_list.take_step();

    forces = Vec3DArray::Zero(p_atoms.n, 3);
    for (auto &pot_ptr : potentials) {
	pot_ptr -> init_step();
	forces += pot_ptr -> calculate_forces(verlet_list);
	pot_ptr -> finish_step(statistics);
    }

    integrator.take_step(p_atoms, forces);
    statistics.take_step();
}

void Simulation::run() {
    while (state.is_running()) {
	take_step();
    }
}

int main(int argc, char **argv) {
    Config config(argc, argv);
    Simulation simulation(config);
    simulation.run();
    return 0;
}

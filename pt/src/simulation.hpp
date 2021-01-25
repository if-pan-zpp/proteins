#pragma once
#include "config.hpp"
#include "state.hpp"
#include "verlist.hpp"
#include "patoms.hpp"
#include "integrator.hpp"
#include "statistics.hpp"
#include "potentials/abstract.hpp"
#include <iostream>
#include <memory>
using namespace std;

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

    Scalar time_delta;
    void calcForces();
};



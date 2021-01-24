#pragma once
#include "def/types.hpp"
#include "conf/constants.hpp"
#include <string>
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
    Config(int argc, char **argv);

    // Angle potential:
    bool angle_pot = true;
    bool dih_angle_pot = false;
    bool tabularized_angle_pot;
    Scalar k_theta = 30.0; // CBA in cg.f

    // Simulation settings:
    int first_index = 0; // first simulated index (nen1 in old code)

    // Temperature constants:
    Scalar delta = 0.005;
    Scalar gamma = 2.0;
    Scalar temp_start = 0.35;
    Scalar temp_end = 0.35;
    Scalar temp_step = 0.05;

    //Predictor-corrector parameters:
    array<Scalar, DER_ORDER + 1>pred_corr_params = {3/16, 251/360, 1, 11/18, 1/6, 1/60};

    // Other configs:
    string out_filename;
    int max_steps = 10;
    Scalar verlet_list_max_eps = 10.0;
};

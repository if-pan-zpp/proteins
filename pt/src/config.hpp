#pragma once
#include "def/types.hpp"
#include "conf/constants.hpp"
#include "conf/units.hpp"
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

    Scalar contact_mltp = 1.3; // MULTIPLIES SIGMA TO CHECK IF CONTACT PRESENT
    bool sink_pot = false; // for sink-like potential for non-native contacts

    // Angle potential:
    bool bond_angle_pot = true;
    bool dih_angle_pot = true;
    bool tabularized_angle_pot = false;
    bool harmonic_dih_pot = false;
    bool pseudo_improper_dihedral_pot = true;
    Scalar k_theta = 30.0;        // CBA in cg.f
    Scalar k_phi_harmonic = 3.33; // CDH in cg.f
    Scalar k_phi_1 = 0.66;        // CDA in cg.f
    Scalar k_phi_2 = 0.66;        // CDB in cg.f

    // PID potential:
    bool mj_matrix = false; // use Miyazawa-Jernigan matrix (lmj in cg.f)
    bool pid_cos = false; // use cos/sin when calculating PID, slower but more precise
    Scalar alpha_bb_pos = 6.4; //alphacos(1,2,3) respectively in cg.f
    Scalar alpha_bb_neg = 6.0;
    Scalar alpha_ss = 1.2;
    Scalar psi0_bb_pos = 1.05;
    Scalar psi0_bb_neg = -1.44;
    Scalar psi0_ss = -0.23;
    Scalar pid_rmin_pos = 5.6_AA;
    Scalar pid_rmin_neg = 6.2_AA;
    Scalar eps_bb = 0.2;
    bool pid_barrier = false;

    // Other potentials:
    Scalar H1 = 50.0_dAA2;
    Scalar H2 = 0.0_dAA4;
    Scalar repulsive_cutoff = 5.0_AA;

    // Simulation settings:
    int first_index = 0; // first simulated index (nen1 in old code)
    Scalar delta = 0.005;

    // Temperature constants:
    Scalar gamma = 2.0;
    Scalar temp_start = 0.35;
    Scalar temp_end = 0.35;
    Scalar temp_step = 0.05;

    // Other configs:
    string out_filename;
    int max_steps = 200;
    Scalar verlet_list_max_eps = 10.0_AA;

    // Cutoffs:
    Scalar all_potentials_r_cut = 18.0_AA;

    // Testing (this is temporary)
    string test_input_file = ""; 
    string test_output_file = "";
};

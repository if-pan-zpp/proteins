#pragma once
#include "def/types.hpp"
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

    // Other configs:
    string out_filename;
    int max_steps = 10;
    def::Scalar verlet_list_max_eps = 10.0;
};
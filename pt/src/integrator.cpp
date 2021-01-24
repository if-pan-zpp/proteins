#include "integrator.hpp"
#include <iostream>
Integrator::Integrator(const Config &_config,
                       const State &_state,
                       const PAtoms &_p_atoms) : state(_state),
                        gamma(_config.gamma), delta(_config.delta),
                        first_index(_config.first_index), residues(_p_atoms.n),
                        pred_corr_params(_config.pred_corr_params) {}

void Integrator::take_step(PAtoms &p_atoms, Vec3DArray &forces) {
    Scalar temperature = state.temperature;

    // lang in cg.f
    const Scalar gamma2 = gamma / delta;
    const Scalar const2 = sqrt(2 * temperature * gamma * delta) * delta;

    const Scalar pseudo_random = 0.28538092181901786; // TODO implement random

    const Scalar pi = acos(-1);

    for(int dim = 0; dim < DIMENSIONS; dim++) {
        for(int res_nr = first_index; res_nr < residues; res_nr++) {
            Scalar r1 = pseudo_random;
            Scalar r2 = pseudo_random;
            Scalar gam = sqrt(-2 * log(r1)) * cos(2 * pi * r2);

            p_atoms.der[1].row(res_nr)[dim] += const2 * gam;

            forces.row(res_nr)[dim] -= gamma2 * p_atoms.der[1](res_nr);
        }
    }
    
    // corr in cg.f
    const Scalar deltsq = 0.5 * delta * delta;

    for(int dim = 0; dim < DIMENSIONS; dim++) {
        for(int res_nr = first_index; res_nr < residues; res_nr++) {
            const Scalar err = p_atoms.der[2].row(res_nr)[dim] - deltsq * forces.row(res_nr)[dim];
            for(int der_nr = 0; der_nr <= DER_ORDER; der_nr++) {
                p_atoms.der[der_nr].row(res_nr)[dim] -= err * pred_corr_params[der_nr];
            }
        }
    }

    //predct in cg.f

    for(int dim = 0; dim < DIMENSIONS; dim++) {
        for(int res_nr = first_index; res_nr < residues; res_nr++) {
            for(int der_nr = 0; der_nr <= DER_ORDER; der_nr++) {
                int mult = 1;
                for(int next_der = der_nr + 1; next_der <= DER_ORDER; next_der++) {
                    mult *= next_der;
                    mult /= max(1, next_der - der_nr);
                    p_atoms.der[der_nr].row(res_nr)[dim] += mult * p_atoms.der[next_der].row(res_nr)[dim];
                }
            }
        }
    }

}

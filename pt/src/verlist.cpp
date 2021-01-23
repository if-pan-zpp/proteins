#include "verlist.hpp"
using namespace std;

VerList::VerList(const Config &_config, const PAtoms &_p_atoms)
    : p_atoms(_p_atoms) {

    eps_upper_bound = _config.verlet_list_max_eps;
}

/*
 * This returns a range in an unmodifiable vector of pairs
 * that contains every pair closer than 'eps'.
 */
pair<VerIt, VerIt> VerList::get_verlet_list(def::Scalar eps) const {
    return {list.begin(), list.end()};
}

/*
 * This calculates the distances between last saved positions and current positions.
 * If it's too big, it recalculates the list.
 */

void VerList::take_step() {
}

/*
 * Potential objects can register the eps that they need before the start of the simulation.
 */
void VerList::register_eps(def::Scalar req_eps) {
    assert (req_eps < eps_upper_bound);
    biggest_req_eps = max(biggest_req_eps, req_eps);
}

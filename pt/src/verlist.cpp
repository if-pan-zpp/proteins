#include "verlist.hpp"
using namespace std;

VerList::VerList(const Config &_config, const PAtoms &_p_atoms)
    : p_atoms(_p_atoms) {

    eps_upper_bound = _config.verlet_list_max_eps;
    reference_pos = Vec3DArray::Zero(_p_atoms.n, 3);
}

/*
 * This returns a range in an unmodifiable vector of pairs
 * that contains every pair closer than 'eps'.
 */
pair<VerIt, VerIt> VerList::get_verlet_list(Scalar eps) const {
    return {list.begin(), list.end()};
}

/*
 * This calculates the distances between last saved positions and current positions.
 * If it's too big, it recalculates the list.
 */

void VerList::take_step() {
    if (need_to_recompute()) {
        list = {};
        Vec3DArray const& pos = p_atoms.der[0];
        auto cutoff_sq = pow(biggest_req_eps, 2);

        // Are klist,kcist, krist etc. updates necessary for the prototype?
        for (int i = 0; i < p_atoms.n; ++i) {
            for (int j = i+1; j < p_atoms.n; ++j) {
                Vec3D dx = pos.row(i) - pos.row(j);
                if (dx.squaredNorm() < cutoff_sq) {
                    list.emplace_back(i, j);
                }
            }
        }
    }
}

bool VerList::need_to_recompute() {
    auto offsets = reference_pos - p_atoms.der[0];
    Eigen::ArrayXd offsets_sqnorm = offsets.matrix().rowwise().squaredNorm();
    Scalar max_allowed_dist = biggest_req_eps / 2; // code_notes.f:660
    auto has_moved_enough = offsets_sqnorm > pow(max_allowed_dist, 2);
    return has_moved_enough.any();
}

/*
 * Potential objects can register the eps that they need before the start of the simulation.
 */
void VerList::register_eps(Scalar req_eps) {
    assert (req_eps < eps_upper_bound);
    biggest_req_eps = max(biggest_req_eps, req_eps);
}

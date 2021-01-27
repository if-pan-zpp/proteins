#include "verlist.hpp"
using namespace std;

VerList::VerList(const Config &_config, const PAtoms &_p_atoms)
    : p_atoms(_p_atoms) {

    eps_upper_bound = _config.verlet_list_max_eps;
    old_positions = Vec3DArray::Zero(_p_atoms.n, 3);

    native_contacts = _p_atoms.native_contacts;
    sort(native_contacts.begin(), native_contacts.end());
    native_contacts.emplace_back(_p_atoms.n + 1, _p_atoms.n + 1);
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
        old_positions = pos;
        vector<FatBool> const& connected = p_atoms.connected;
        Scalar cutoff_sq = pow(biggest_req_eps + verlet_cutoff, 2);

        int native_it = 0;
        for (int i = 0; i < p_atoms.n; ++i) {
            int start_j = i + 2;
            if (connected[i] && connected[i + 1]) start_j++;
            for (int j = start_j; j < p_atoms.n; ++j) {
                Vec3D dx = pos.row(i) - pos.row(j);
                if (dx.squaredNorm() < cutoff_sq) {
                    const pair<int, int> p{i, j};
                    while (native_contacts[native_it] < p) {
                        native_it++;
                    }

                    if (native_contacts[native_it] != p) {
                        list.push_back(p);
                    }
                }
            }
        }
    }
}

bool VerList::need_to_recompute() {
    auto offsets = old_positions - p_atoms.der[0];
    Eigen::ArrayXd offsets_sqnorm = offsets.matrix().rowwise().squaredNorm();
    Scalar max_allowed_dist = verlet_cutoff / 2; // code_notes.f:660
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

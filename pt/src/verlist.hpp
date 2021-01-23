#pragma once
#include "config.hpp"
#include "patoms.hpp"
#include "def/types.hpp"
using namespace std;

using VerIt = vector<pair<int,int>>::const_iterator;

/*
 * Class that provides list of pairs of pseudoatoms that are close.
 * The exact design and implementation is a very hard and interesting problem.
 * But for the purposes of the prototype that follows cg.f in its implementation,
 * it's easy for now.
 */
class VerList {
public:
    VerList(const Config &_config, const PAtoms &_p_atoms);

    /*
     * This returns a range in an unmodifiable vector of pairs
     * that contains every pair closer than 'eps'.
     */
    pair<VerIt, VerIt> get_verlet_list(def::Scalar eps) const;

    /*
     * This calculates the distances between last saved positions and current positions.
     * If it's too big, it recalculates the list.
     */
    void take_step();

    /*
     * Potential objects can register the eps that they need before the start of the simulation.
     */
    void register_eps(def::Scalar req_eps);

private:
    const PAtoms &p_atoms;
    def::Scalar eps_upper_bound;
    def::Scalar biggest_req_eps = 0.;
    vector<pair<int, int>> list;
};

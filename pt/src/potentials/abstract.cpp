/*
 * This file defines the list of all supported potentials.
 */
#include "abstract.hpp"

#include "bondAngle.hpp"
#include "dihAngle.hpp"
#include "structuredLJ.hpp"
#include "localRepulsive.hpp"
#include "harmonicTethering.hpp"
#include "globalRepulsive.hpp"
#include "pseudoImproperDihedral.hpp"

vector<unique_ptr<Potential>> all_supported_potentials(const Config &config,
                                                       const State &state,
                                                       const PAtoms &p_atoms) {
    vector<unique_ptr<Potential>> all;
    all.emplace_back(make_unique<BondAngle>(config, state, p_atoms));
    all.emplace_back(make_unique<DihAngle>(config, state, p_atoms));
    all.emplace_back(make_unique<StructuredLJ>(config, state, p_atoms));
    all.emplace_back(make_unique<LocalRepulsive>(config, state, p_atoms));
    all.emplace_back(make_unique<HarmonicTethering>(config, state, p_atoms));
    all.emplace_back(make_unique<GlobalRepulsive>(config, state, p_atoms));
    all.emplace_back(make_unique<PseudoImproperDihedral>(config, state, p_atoms));
    return all;
}

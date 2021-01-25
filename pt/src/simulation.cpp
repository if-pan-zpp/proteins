#include "simulation.hpp"

Simulation::Simulation(const Config &config) :
    state(config),
    p_atoms(config),
    verlet_list(config, p_atoms),
    integrator(config, state, p_atoms),
    statistics(config, state, p_atoms) {

    for (auto &pot_ptr: all_supported_potentials(config, state, p_atoms)) {
        if (pot_ptr -> is_enabled()) potentials.emplace_back(move(pot_ptr));
    }

    for (auto &pot_ptr : potentials) pot_ptr -> init_and_register(verlet_list);

    forces = Vec3DArray::Zero(p_atoms.n, 3);
    
    time_delta = config.delta;
}

void Simulation::calcForces() {
    verlet_list.take_step();

    forces = Vec3DArray::Zero(p_atoms.n, 3);
    for (auto &pot_ptr : potentials) {
        pot_ptr -> init_step();
        forces += pot_ptr -> calculate_forces(verlet_list);
        pot_ptr -> finish_step(statistics);
    }
}

void Simulation::run() {
    calcForces();
    const Scalar time_delta_sq = 0.5 * time_delta * time_delta;
    p_atoms.der[2] = time_delta_sq * forces;

    while (state.is_running()) {
        state.take_step();
        integrator.take_step(p_atoms, forces);

        calcForces();
    
        statistics.take_step();
    }
    cout << "POSITIONS:\n";
    cout << p_atoms.der[0] << endl;
}

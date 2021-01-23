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

    cout << "Number of enabled potentials: " << potentials.size() << endl;
}

void Simulation::take_step() {
    state.take_step();
    verlet_list.take_step();

    forces = Vec3DArray::Zero(p_atoms.n, 3);
    for (auto &pot_ptr : potentials) {
        pot_ptr -> init_step();
        forces += pot_ptr -> calculate_forces(verlet_list);
        pot_ptr -> finish_step(statistics);
    }

    integrator.take_step(p_atoms, forces);
    statistics.take_step();
}

void Simulation::run() {
    while (state.is_running()) {
        take_step();
    }
}

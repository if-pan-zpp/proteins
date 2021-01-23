#include "config.hpp"
#include "simulation.hpp"
using namespace std;

int main(int argc, char **argv) {
    Config config(argc, argv);
    Simulation simulation(config);
    simulation.run();
    return 0;
}

#include "config.hpp"

Config::Config(int argc, char **argv) {
    // Here there will be a lot of parsing configs, checking them,
    // printing to the user, prompting the user...

    if (argc > 1) {
        test_input_file = string(argv[1]);
    }
    if (argc > 2) {
        test_output_file = string(argv[2]);
    }
}

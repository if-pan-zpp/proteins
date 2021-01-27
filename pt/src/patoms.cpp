#include "patoms.hpp"

#include <iostream>
#include <fstream>

PAtoms::PAtoms(const Config &config) {
    n = 100;
    der[0] = Eigen::ArrayX3d::Random(n, 3);
    for (int i = 1; i<= DER_ORDER; i++) {
        der[i] = Eigen::ArrayX3d::Zero(n, 3);
    }
    connected.resize(n, 1);
    native_pos = Eigen::ArrayX3d::Random(n, 3);

    // Read initial conditions from a test file.
    if (!config.test_input_file.empty()) {
        n = -1;
        ifstream inp(config.test_input_file, ifstream::in);

        string word;
        Scalar x, y, z;

        while (inp >> word) {
            if (n == -1) assert (word == "n");

            if (word == "n") {
                inp >> n;
                for (int i = 0; i <= DER_ORDER; ++i)
                    der[i] = Eigen::ArrayX3d::Zero(n, 3);
                connected.resize(n, 1);
                native_pos = Eigen::ArrayX3d::Zero(n, 3);
            }
            else if (word == "rng_calls_cnt") {
                inp >> rng_calls_cnt;
            }
            else if (word == "native_pos") {
                for (int i = 0; i < n; ++i) {
                    inp >> x >> y >> z;
                    native_pos.row(i) = Eigen::Array3d(x, y, z);
                }
            }
            else if (word == "positions") {
                for (int i = 0; i < n; ++i) {
                    inp >> x >> y >> z;
                    der[0].row(i) = Eigen::Array3d(x, y, z);
                }
            }
            else if (word == "velocities") {
                for (int i = 0; i < n; ++i) {
                    inp >> x >> y >> z;
                    der[1].row(i) = Eigen::Array3d(x, y, z);
                }
            }
            else if (word == "native_contacts") {
                int cnt;
                inp >> cnt;
                for (int i = 0; i < cnt; ++i) {
                    int a, b;
                    inp >> a >> b;
                    native_contacts.emplace_back(a - 1, b - 1);
                }
            }
        }
    }
}

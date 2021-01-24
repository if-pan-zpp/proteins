#include "patoms.hpp"

PAtoms::PAtoms(const Config &config) {
    n = 100;
    der[0] = Eigen::ArrayX3d::Random(n, 3);
    for (int i = 1; i<= DER_ORDER; i++)
        der[i] = Eigen::ArrayX3d::Zero(n, 3);
    connected.resize(n, 1);
    native_theta.resize(n, 0.5);
}

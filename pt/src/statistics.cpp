#include "statistics.hpp"
#include <iostream>
#include <fstream>

Statistics::Statistics(const Config &_config,
                       const State &_state,
                       const PAtoms &_p_atoms)
    : state(_state),
    p_atoms(_p_atoms),
    test_output_file(_config.test_output_file) {};

void Statistics::take_step() {
    cerr << "We're in step " << state.cur_step << endl;
}

void Statistics::print_summary() {
    cerr << "POSITIONS:\n";
    cerr << p_atoms.der[0] << endl;

    // JUST FOR TESTING:
    if (!test_output_file.empty()) {
        ifstream inp(test_output_file, ifstream::in);

        string field;
        inp >> field;
        assert (field == "n");
        int n;
        inp >> n;
        assert (n == p_atoms.n);

        Vec3DArray ref_positions = Vec3DArray::Zero(n, 3);
        for (int i = 0; i < n; ++i) {
            double x, y, z;
            inp >> x >> y >> z;
            ref_positions.row(i) = Vec3D(x, y, z);
        }

        Vec3DArray pos_diff = p_atoms.der[0] - ref_positions;
        Eigen::ArrayXd diff_lengths = pos_diff.matrix().rowwise().norm();
        cout << "Max diff = " << diff_lengths.maxCoeff() << endl;
        cout << "Average diff = " << diff_lengths.mean() << endl;
    }
}

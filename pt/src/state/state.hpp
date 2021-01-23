#pragma once
#include "def.hpp"
#include "ir.pb.hpp"
#include <string>

namespace state {
    class State {
    public:
        VectorList pos, vel;
        ScalarList mass;
        std::string seq;

        static State fromIR(pdb::ir::Entry const& entry);

        void force(VectorList& f);
    };
}

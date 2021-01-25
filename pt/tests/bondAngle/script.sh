#!/bin/bash

# The number of steps ./cg simulates is 200*mstep (set in inputfile1)
# By default, ./pt makes 200 steps, this can be changes in Config class now.

REF_PATH=../../../ref/

echo "Before using, copy executable pt here."
cp $REF_PATH/code/cg.f .
cp $REF_PATH/data/example1/1ubq.pdb .
patch cg.f cg_patch -o cg_mod.f
gfortran -O3 cg_mod.f -o cg
./cg inputfile1
./pt test_input.txt
echo "Compare with test1.pdb (for TIME=mstep, not TIME=0)"

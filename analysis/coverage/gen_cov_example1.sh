#!/bin/bash

cp ../f77-source/code/cg.f .
cp ../f77-source/data/example1/* .
# Decrease numbers so that it finishes quickly
sed 's/ntraj 10/ntraj 1/' inputfile1 | sed 's/mstep .*/mstep 150/' > inputfile
rm inputfile1
make

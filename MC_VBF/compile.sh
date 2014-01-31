#! /bin/bash

# standard "rivet-buildplugin Rivet_VBF.so MC_VBF.cc" only with extra "-L/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67/rivet/2.0.0/x86_64-slc6-gcc47-opt/lib"

#/afs/cern.ch/sw/lcg/contrib/gcc/4.7.2/x86_64-slc6-gcc47-opt/bin/c++ -o "Rivet_VBF.so" -shared -fPIC -L/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67/rivet/2.0.0/x86_64-slc6-gcc47-opt/lib -I/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67/rivet/2.0.0/x86_64-slc6-gcc47-opt/include -I/afs/cern.ch/sw/lcg/external/HepMC/2.06.08/x86_64-slc6-gcc47-opt/include -I/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67/yoda/1.0.4/x86_64-slc6-gcc47-opt/include -I/afs/cern.ch/sw/lcg/external/fastjet/3.0.3/x86_64-slc6-gcc47-opt/include -I/afs/cern.ch/sw/lcg/external/GSL/1.10/x86_64-slc6-gcc47-opt/include -I/afs/cern.ch/sw/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc47-opt/include/boost-1_53 -std=c++11 -pedantic -ansi -Wall -Wno-long-long -Wno-format -Werror=uninitialized -Werror=delete-non-virtual-dtor -O2 -Wl,--no-as-needed -L/lib -L/afs/cern.ch/sw/lcg/external/HepMC/2.06.08/x86_64-slc6-gcc47-opt/lib -L/afs/cern.ch/sw/lcg/external/MCGenerators_lcgcmt67/yoda/1.0.4/x86_64-slc6-gcc47-opt/lib -Wl,-rpath,/afs/.cern.ch/sw/lcg/external/fastjet/3.0.3/x86_64-slc6-gcc47-opt/lib -lm -L/afs/.cern.ch/sw/lcg/external/fastjet/3.0.3/x86_64-slc6-gcc47-opt/lib -lfastjettools -lfastjet -lfastjetplugins -lsiscone_spherical -lsiscone -L/afs/.cern.ch/sw/lcg/contrib/gcc/4.7.2/x86_64-slc6-gcc47-opt/bin/../lib/gcc/x86_64-unknown-linux-gnu/4.7.2 -L/afs/.cern.ch/sw/lcg/contrib/gcc/4.7.2/x86_64-slc6-gcc47-opt/bin/../lib/gcc -L/afs/.cern.ch/sw/lcg/contrib/gcc/4.7.2/x86_64-slc6-gcc47-opt/bin/../lib/gcc/x86_64-unknown-linux-gnu/4.7.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/afs/.cern.ch/sw/lcg/contrib/gcc/4.7.2/x86_64-slc6-gcc47-opt/bin/../lib/gcc/x86_64-unknown-linux-gnu/4.7.2/../../.. -lgfortran -lm -lquadmath MC_VBF.cc -lRivet
rivet-buildplugin Rivet_VBF.so MC_VBF.cc
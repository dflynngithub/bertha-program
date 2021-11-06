#!/usr/bin/env bash
NAME=${1?Error: no name given}
#ulimit -s unlimited 
#./bertha-budget < input/$NAME
#./bertha-nuclear < input/$NAME
# gfortran -O4 bertha-atomic.f -o bertha-atomic -llapack -lblas
#./bertha-atomic < input/$NAME
gfortran -O4 bertha-scf.f -o bertha-scf -llapack -lblas
./bertha-scf < input/$NAME
#./bertha-mbpt < input/$NAME
#./bertha-rspt < input/$NAME
#./bertha-visual < input/$NAME
#
# To run this with an input file 'input':
# > ./run-bertha.sh input
#
# In the following:
# If a space is required, it will be written explicitly in the option.
# Options with multiple choices are denoted with [ ], but don't actually use brackets.
#
# Compiler options before 'gfortran'
# ulimit -a            : see current memory usage limits
# ulimit -[ ] unlimited: for dealing with high memory usage, with some cases:
#          d           : data segment size (kB)
#          f           : file sizes (blocks)
#          l           : maximum locked memory (kB)
#          m           : maximum memory size (kB)
#          n           : open files (minimum value)
#          s           : stack size (kB)
#          t           : CPU time (seconds)
#          u           : maximum user processes
#          v           : virtual memory (kB)
#
# Compiler options between 'gfortran' and 'bertha.f'
# -mcmodel=[      ]: extend static memory allocation if COMMON arrays have size over 2GB
#           small  : restrict code and data to first 2GB of address space (default)
#           medium : restrict compiler to the first 2GB but place no memory restriction on data
#           large  : place no memory restriction on code or data
# -O[ ]            : optimisation level N={1,2,3,4,5} (compiler takes longer but is worthwhile)
#    1             : minimum of statement-level optimisations
#    2             : basic block level optimisations (gives the smallest code size)
#    3             : adds loop unrolling and global optimisations at the function level
#    4             : adds automatic inlining of routines contained in the same file
#    5             : attempt aggressive optimisations (might be dangerous)
# -g               : compile for debugging
# -time            : time each compilation phase (look at importing result to BERTHA)
#
# Compiler options after 'bertha.f'
# -o [filename]: specify executable name
# -l[x]        : add library 'libx.a' to linker's list of search libraries
#    lapack    : linear algebra library
#    blas      : basic linear algebra sub-programs
#    fftw3     : Fourier transform software
#    mpi       : OpenMPI library
# -fopenmp     : link the OpenMP parallelisation option
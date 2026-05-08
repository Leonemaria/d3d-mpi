#!/bin/bash
export OMP_NUM_THREADS=$2
mpiexec -n $1 -bind-to numa ./d3d

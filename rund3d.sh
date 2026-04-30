#!/bin/bash
export OMP_NUM_THREADS=$2
export OMP_PLACES=cores
export OMP_PROC_BIND=close
mpiexec -n $1 -bind-to numa ./d3d

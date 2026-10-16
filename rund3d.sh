#!/bin/bash
mpiexec -n $1 --bind-to numa --map-by numa ./d3d $2

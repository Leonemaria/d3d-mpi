#!/bin/bash
mpiexec -n $1 -bind-to numa ./d3d $2

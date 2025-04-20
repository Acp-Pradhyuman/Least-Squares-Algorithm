#!/bin/bash

# Name of the compiled MPI executable
# EXEC=./mpi_least_squares
EXEC=./a.out

# Check if executable exists
if [ ! -f "$EXEC" ]; then
    echo "Executable $EXEC not found!"
    exit 1
fi

# Loop through powers of 2 from 1 to 64
for ((i=0; i<=6; i++)); do
    NP=$((2**i))
    echo "Running with $NP MPI processes..."
    mpirun -np $NP $EXEC
    echo ""
done
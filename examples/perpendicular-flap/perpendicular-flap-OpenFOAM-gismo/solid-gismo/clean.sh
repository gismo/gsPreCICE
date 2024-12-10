#!/bin/bash

# Cleaning script for solid-gismo-elasticity directory

echo "Cleaning unnecessary files..."

# Remove precice-profiling directory
if [ -d "precice-profiling" ]; then
    rm -rf precice-profiling
    echo "Deleted 'precice-profiling' folder."
fi

# Remove precice-run directory
if [ -d "../precice-run" ]; then
    rm -rf ../precice-run
    echo "Deleted 'precice-run' folder."
fi

# Remove files ending with .pvd, .vts, .vtp, .log, and .txt
for ext in pvd vts vtp log txt; do
    find . -type f -name "*.$ext" -exec rm -f {} \;
    echo "Deleted all *.$ext files."
done

echo "Cleaning completed!"

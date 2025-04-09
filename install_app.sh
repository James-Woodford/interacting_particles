#!/bin/bash
set -e
echo "Configuring and building the project..."

# Create build directory.
if [ -d "build" ]; then
    echo "Removing old build directory..."
    rm -rf build
fi
mkdir build
cd build

# Run CMake
echo "Running CMake configuration..."

PYTHON_EXEC=$(which python3.10)
cmake .. -DPYTHON_EXECUTABLE=$PYTHON_EXEC

# Build
echo "Building the project..."
cmake --build .
cd ..
cp build/interacting_particles.cpython-310-x86_64-linux-gnu.so .
echo "Build complete"

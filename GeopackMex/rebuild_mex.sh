#!/bin/bash

# ==============================================================================
# Build Script for IGRF MEX File
# Description: Compiles Fortran (Geopack/TS04) and C sources into a MATLAB MEX file.
# Platform: macOS (Apple Silicon / ARM64)
# Dependencies: MATLAB, gfortran (via Homebrew)
# ==============================================================================

# 1. Clean up previous build artifacts and temporary directories
rm -rf tmp/
rm igrfmex.mexmaca64

# 2. Configuration: Set Paths
# ------------------------------------------------------------------------------
# Define MATLAB installation root path.
# NOTE: Update this path according to your local MATLAB installation.
MATLABROOT="/Applications/MATLAB_R2024b.app"

# Locate the gfortran compiler (Assumes installation via Homebrew).
GFORT=$(which gfortran)

# Define the path to the MATLAB 'mex' command executable.
MEX="$MATLABROOT/bin/mex"

# Retrieve the installation prefix for GCC/gfortran to locate libraries.
GCC_PREFIX=$(brew --prefix gcc)

# 3. Preparation: Create Temporary Directory
# ------------------------------------------------------------------------------
# Create a directory to store intermediate object files (.o).
TMP_DIR="tmp"
mkdir -p "$TMP_DIR"

# 4. Compilation: Fortran Source Files
# ------------------------------------------------------------------------------
# Compile Geopack and TS04 Fortran codes into object files.
echo "Compiling Fortran sources..."
MACOSX_DEPLOYMENT_TARGET=11.0 $GFORT -O3 -arch arm64 -o "$TMP_DIR/geopack.o" -c "./Geopack-2008_dp.for"
MACOSX_DEPLOYMENT_TARGET=11.0 $GFORT -O3 -arch arm64 -o "$TMP_DIR/ts04.o" -c "./TS04c.for"

# 5. Compilation: C Interface
# ------------------------------------------------------------------------------
# Compile the C wrapper script using MATLAB's mex command.
# NOTE: The library path 'lib/gcc/14' depends on the installed GCC version (e.g., 12, 13, 14).
echo "Compiling C interface..."
$MEX -R2018a -O -outdir "$TMP_DIR" -c "./igrfmex.c" -L"$GCC_PREFIX/lib/gcc/14" -L"/usr/local/lib" -lgfortran -lquadmath

# 6. Linking: Generate Final MEX Binary
# ------------------------------------------------------------------------------
# Link all object files and libraries to generate the final .mexmaca64 file.
echo "Linking and generating MEX file..."
$MEX -R2018a -outdir "." "$TMP_DIR/igrfmex.o" "$TMP_DIR/geopack.o" "$TMP_DIR/ts04.o"  -L"$GCC_PREFIX/lib/gcc/14" -L"/usr/local/lib" -lgfortran -lquadmath

echo "Build complete."
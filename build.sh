#!/bin/bash
# Download LAMMPS, add the RBC styles of src/ and build lmp for one
# NVIDIA GPU (KOKKOS package with CUDA, no MPI).
#
#   ./build.sh                  # GPU type from nvidia-smi
#   ARCH=AMPERE86 ./build.sh    # or name it (Kokkos architecture):
#                               #   VOLTA70   V100
#                               #   TURING75  RTX 2080, T4
#                               #   AMPERE80  A100
#                               #   AMPERE86  RTX 3090, A6000, A40
#                               #   ADA89     L40S, RTX 4090, RTX 6000 Ada
#                               #   HOPPER90  H100, H200
#   NPROC=16 ./build.sh         # more parallel compile jobs (default 4)
#
# Needs git, CMake 3.20 or newer, a C++17 compiler and the CUDA toolkit.
# On a cluster, load them first, for example
#   module load gcc cuda cmake
set -e

HERE=$(cd "$(dirname "$0")" && pwd)
TAG=${LAMMPS_TAG:-stable_30Sep2026}
LAMMPS=${LAMMPS_DIR:-$HERE/lammps}
NPROC=${NPROC:-4}

if [ -z "$ARCH" ]; then
  cap=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d ' .')
  case "$cap" in
    70) ARCH=VOLTA70 ;;
    75) ARCH=TURING75 ;;
    80) ARCH=AMPERE80 ;;
    86) ARCH=AMPERE86 ;;
    89) ARCH=ADA89 ;;
    90) ARCH=HOPPER90 ;;
    *) echo "No GPU found here, or an unknown one: set ARCH (see the top of build.sh)."; exit 1 ;;
  esac
fi
echo "building for Kokkos architecture $ARCH"

if [ ! -d "$LAMMPS" ]; then
  git clone --depth 1 --branch "$TAG" https://github.com/lammps/lammps.git "$LAMMPS"
fi
cp "$HERE"/src/*.cpp "$HERE"/src/*.h "$LAMMPS"/src/
cp "$HERE"/src/KOKKOS/*.cpp "$HERE"/src/KOKKOS/*.h "$LAMMPS"/src/KOKKOS/

cmake -S "$LAMMPS"/cmake -B "$HERE"/build -C "$LAMMPS"/cmake/presets/kokkos-cuda.cmake \
  -D CMAKE_BUILD_TYPE=Release -D Kokkos_ARCH_"$ARCH"=yes -D BUILD_MPI=no -D BUILD_OMP=no \
  -D PKG_MOLECULE=yes -D PKG_DPD-BASIC=yes -D PKG_EXTRA-FIX=yes
cmake --build "$HERE"/build -j "$NPROC"
echo "built $HERE/build/lmp"

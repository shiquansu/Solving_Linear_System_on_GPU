#!/bin/bash
echo "reset module environment"
source /etc/profile.d/modules.sh
module purge
module load intel/oneapi/2023.2.0
module load intel-comp-rt/agama-ci-devel
#module load intel/pti-gpu/2022-07-06 
module load intel/pti-gpu-nda
module load intel/vtune/2023.2.0

echo "MKLROOT="${MKLROOT}
vtune --version

echo "LD_LIBRARY_PATH: "${LD_LIBRARY_PATH}

echo "CPATH: "$CPATH

echo "clean up *.mod *.modmic"
rm *.mod *.modmic core
echo "compile:"

export LIBOMPTARGET_PLUGIN_PROFILE=T
make clean
make
#ifx -DDC -i8 -DMKL_ILP64 -qopenmp -fopenmp-targets=spir64 -fsycl -free solve_complex_matrix.f90 -o solve_complex_matrix -L${MKLROOT}/lib/intel64 -lmkl_sycl -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl

echo "run:"
#OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0.1 LIBOMPTARGET_PLUGIN_PROFILE=T LIBOMPTARGET_DEBUG=1 ./solve_complex_matrix -n 1424 -b 2 -r 1 -c 2 >& solve_complex_matrix.out 


MKL_VERBOSE=2 OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0.1 ./solve_complex_matrix

MKL_VERBOSE=2 OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0.1 vtune -c gpu-offload ./solve_complex_matrix

MKL_VERBOSE=2 OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0.1 vtune -c gpu-hotspots ./solve_complex_matrix

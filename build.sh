#!/bin/bash
echo "reset module environment"
source /etc/profile.d/modules.sh
module purge
module load intel/oneapi/2023.2.0 # for openmp dispatch approach
#module load intel/oneapi/2023.1.0 # for c wrapper approach
module load intel-comp-rt/agama-ci-devel
#module load intel/pti-gpu/2022-07-06 
#module load intel/pti-gpu-nda
#module load intel/vtune/2023.2.0

echo "MKLROOT="${MKLROOT}
vtune --version

#echo "LD_LIBRARY_PATH: "${LD_LIBRARY_PATH}
#echo "CPATH: "$CPATH

echo "clean up *.mod *.modmic"
rm *.mod *.modmic core
echo "compile:"

#mkl dispatch performance related flags:
#https://www.intel.com/content/www/us/en/docs/oneapi/optimization-guide-gpu/2023-2/prefetching.html
export IGC_allowLICM=0
unset IGC_DisableMatchMad
export EnableImplicitScaling=1
export IGC_ForceOCLSIMDWidth=16



export LIBOMPTARGET_PLUGIN_PROFILE=T
LIBOMPTARGET_DEBUG=0

make clean
make scm

echo "run:"

MKL_VERBOSE=2 OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0.1 ./scm

#OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0.1 vtune -c gpu-offload ./scm

#IGC_allowLICM=0 IGC_ForceOCLSIMDWidth=32 OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0.1 ./scm

#IGC_allowLICM=0 IGC_ForceOCLSIMDWidth=32 OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0.1 vtune -c gpu-hotspots ./scm 2>&1 > screenoutput_gh.txt
#For onetrace and vtune on how many tile/slice onemkl dispatch calls are using: I recommend talking to Ma, Zhiqiang <zhiqiang.ma@intel.com> and Caday, Peter <peter.caday@intel.com>
# https://www.intel.com/content/www/us/en/docs/oneapi/optimization-guide-gpu/2023-1/intel-profiling-tools-interfaces-for-gpu.html
#OMP_TARGET_OFFLOAD=MANDATORY ZE_AFFINITY_MASK=0.1 onetrace --chrome-call-logging --chrome-device-timeline ./scm 2>&1 > screenoutput_dt.txt

#! /bin/bash
#
# Submit job as (build defaults to Release):
#
#   sbatch compile.sh
#   sbatch --export=BUILD='Debug',ALL compile.sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=2:00:00
#SBATCH --mem=32GB

# Load modules
module purge
module load mercurial/intel/4.0.1
module load gcc/6.3.0
module load cmake/intel/3.7.1
module load mpfr/gnu/3.1.5
module load gmp/gnu/6.1.2

export GMP_DIR=$GMP_ROOT
export MPFR_DIR=$MPFR_ROOT

export CC=${GCC_ROOT}/bin/gcc
export CXX=${GCC_ROOT}/bin/g++

# Microstructure code
module load suitesparse/intel/4.5.4
module load boost/intel/1.62.0
module load glog/intel/0.3.4

# Run job
cd "${SLURM_SUBMIT_DIR}"

# Compile main program
SOURCE_DIR=${SLURM_SUBMIT_DIR}
#BUILD_DIR=build
BUILD_DIR=/scratch/${USER}/build/microstructures

mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

echo ${BUILD}
if [ -z "${BUILD}" ]; then
	BUILD=Release
fi

mkdir -p ${BUILD}
pushd ${BUILD}
cmake -DCMAKE_BUILD_TYPE=${BUILD} -DMESHFEM_BUILD_BINARIES=ON -DCMAKE_CXX_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=0" ${SOURCE_DIR}
make -j20
make test
popd

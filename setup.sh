#!/bin/bash
source /cvmfs/fcc.cern.ch/sw/latest/setup.sh
#source /cvmfs/sw.hsf.org/spackages/linux-centos7-broadwell/gcc-8.3.0/key4hep-stack-2021-03-26-ch6gml3x3wk67ydbiekh7yx7bmr7nmil/setup.sh
#source /cvmfs/sw.hsf.org/contrib/spack/share/spack/setup-env.sh

export PYTHONPATH=$PWD:$PYTHONPATH
export LD_LIBRARY_PATH=$PWD/install/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$PWD/install/include/FCCAnalyses:$ROOT_INCLUDE_PATH

#spack load acts@5.0.0 arch=linux-centos7-x86_64

#spack load acts@5.00.0 build_type=Debug
#export LOCALFCCANALYSES=$PWD/install/include/FCCAnalyses
#spack load py-pyyaml

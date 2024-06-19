#!/usr/bin/env bash

# Script to setup key4hep nightly envrionment that allows to use (py)torch, onnx and Marlin
# together with DD4hep.
#
#


source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
# source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh -r 2024-03-07

TORCH_PATH=$(dirname $(python -c 'import torch; print(f"{torch.__file__}")'))
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}:${TORCH_PATH}/share/cmake
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${TORCH_PATH}/lib

ONNXRUNTIME_PATH=$(dirname $(python -c 'import onnxruntime; print(f"{onnxruntime.__file__}")'))
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(realpath ${ONNXRUNTIME_PATH}/../../../../)/lib64  # CentOS7
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(realpath ${ONNXRUNTIME_PATH}/../../../../)/lib    # Ubuntu

#!/bin/bash

source ./miniconda3/bin/activate rascaline_install

export CONDA_BASE="$PWD/miniconda3/bin"

export OMP_NUM_THREADS=1

ipi=i-pi



plumed_path=$(which plumed)
echo "Plumed is located at: $plumed_path"

i-pi input.xml &> log.i-pi & 
sleep 2
python run-mace.py -u -a driver-works-100 -m driver-works-100
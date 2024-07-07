#!/bin/bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/eos/user/w/wochung/src/SCEPCAL/install/lib64
export PYTHONPATH=$PYTHONPATH:/eos/user/w/wochung/src/SCEPCAL/install/python

ddsim /eos/user/w/wochung/src/SCEPCAL/scripts/scepcal_steering.py
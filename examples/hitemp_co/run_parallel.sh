#!/bin/bash

CONFIG_FILE=hitemp_co.py

mkdir logs
# Run 4 pyROX-calls in parallel with different temperatures
pyROX $CONFIG_FILE -c -pbar --T_grid 500 -out xsec_T500K.hdf5 > logs/range_T500K.log 2>&1 &
pyROX $CONFIG_FILE -c -pbar --T_grid 1000 -out xsec_T1000K.hdf5 > logs/range_T1000K.log 2>&1 &
pyROX $CONFIG_FILE -c -pbar --T_grid 2000 -out xsec_T2000K.hdf5 > logs/range_T2000K.log 2>&1 &
pyROX $CONFIG_FILE -c -pbar --T_grid 3000 -out xsec_T3000K.hdf5 > logs/range_T3000K.log 2>&1 &
wait
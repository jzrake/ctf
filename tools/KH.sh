#!/bin/bash

resolutions=(32 64 128 256)
scheme=hllc-plm-rk3
mpi="mpiexec -np 2"

for res in ${resolutions[@]}
do
    echo "******************* Running $res^2 *******************"
    $mpi ./ctf examples/KH.lua \
	--tmax=2.5 \
	--resolution=$res \
	--advance=rk3 \
	--riemann=hllc \
	--solver=godunov \
	--id=$res-$scheme \
	--cpi=0.025 \
	--output=KH-$res-$scheme.h5
done

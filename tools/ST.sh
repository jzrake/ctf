#!/bin/bash

problems=(
    'Shocktube1'
    'Shocktube2'
    'Shocktube3'
    'Shocktube4'
    'Shocktube5'
    'ContactWave')

for p in ${problems[@]}
do
    echo "******************* Running $p *******************"
    ./ctf examples/tests-1d.lua $p \
	--resolution=128 \
	--id=hllc-plm-rk3 \
	--backend=mara \
	--output=$p.h5
    ./tools/plot1d.py $p.h5 -o $p.pdf
    ./tools/plot1d.py $p.h5 -o $p.png
done

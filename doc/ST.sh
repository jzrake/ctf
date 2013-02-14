

declare -a problems=('Shocktube1' 'Shocktube2' 'Shocktube3' 'Shocktube4' 'Shocktube5')

for i in {0..4}
do
    echo ${problems[i]}
    ./ctf examples/tests-1d.lua ${problems[i]} \
	--id=hllc-plm-rk3 \
	--backend=mara \
	--output=${problems[i]}.h5
done


# Give the problem name as the first argument
problem=$1

# ----------------------------------------------------------
# Mara: HLLC-PLM-MUSCL
# ----------------------------------------------------------
./ctf examples/tests-1d.lua ${problem} --convergence \
    --backend=mara \
    --tmax=0.5 \
    --solver=muscl \
    --advance=single \
    --output=hllc-plm-muscl.dat

# ----------------------------------------------------------
# Mara: HLLC-PLM-RK3
# ----------------------------------------------------------
./ctf examples/tests-1d.lua ${problem} --convergence \
    --backend=mara \
    --tmax=0.5 \
    --solver=godunov \
    --advance=rk3 \
    --output=hllc-plm-rk3.dat

# ----------------------------------------------------------
# Fish: HLLC-WENO5-RK3
# ----------------------------------------------------------
./ctf examples/tests-1d.lua ${problem} --convergence \
    --backend=fish \
    --tmax=0.5 \
    --reconstruction=weno5 \
    --solver=godunov \
    --advance=rk3 \
    --output=hllc-weno5-rk3.dat

# ----------------------------------------------------------
# Mara: CHAR-WENO5-RK3
# ----------------------------------------------------------
./ctf examples/tests-1d.lua ${problem} --convergence \
    --backend=mara \
    --tmax=0.5 \
    --solver=spectral \
    --advance=rk3 \
    --output=char-weno5-rk3.dat

# ----------------------------------------------------------
# Mara: CHAR-WENO5-RK4
# ----------------------------------------------------------
./ctf examples/tests-1d.lua ${problem} --convergence \
    --backend=mara \
    --tmax=0.5 \
    --solver=spectral \
    --advance=rk4 \
    --output=char-weno5-rk4.dat


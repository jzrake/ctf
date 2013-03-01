#!/bin/bash

if [ "$1" == "--explain" ]; then
    echo "Run a suite of self-tests"
    exit 0
fi

ctf=bin/ctf-main

$ctf --splash | grep HDF5 | grep -q yes && $ctf lua-hdf5/run.lua
$ctf --splash | grep HDF5 | grep -q yes && $ctf lua-hdf5/LuaHDF5.lua
$ctf --splash | grep MPI  | grep -q yes && $ctf lua-mpi/run.lua

$ctf fish/smr1d.lua

$ctf modules/class.lua
$ctf modules/array.lua
$ctf modules/mesh.lua
$ctf modules/FishClasses.lua

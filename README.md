
# Overview

CTF stands for Computational Turbulence Framework. That may be a little
ambitious. Really it's a just hydrodynamics code with lots of room for
expansion. It currently includes the following projects as submodules, each of
them exposed as a Lua module:

+ `Mara`: Relativistic MHD turbulence code
+ `C.O.W.`: 'Cube of Wonder' analysis tools
+ `lua-hdf5`: High-level parallel I/O with HDF5
+ `lua-mpi`: Lua interface to the Message Passing Interface
+ `lua-glut`: Lua wrappers for OpenGL and GLUT


# Build instructions

Copy the `Makefile.in.template` file to `Makefile.in`, and modify it for your
system.

+ `git submodule update --init`
+ `make lua`
+ `make`


# License

This code is made freely available for anybody's use. I only ask that if you
find it useful, please send me an email about how it works, and in what project
you are using it.


The CTF is licensed under the same terms as Lua itself.

Copyright (c) 2012, Jonathan Zrake <jonathan.zrake@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


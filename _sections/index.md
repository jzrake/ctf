

### Overview
#+{'figures':
#+[{'src': 'turbulent-bndlayer-KH-mara-small.png',
#+  'url': 'http://zrake.webfactional.com/static/images/turbulent-bndlayer-KH-mara.png',
#+  'alt': 'Kelvin-Helmholtz instability: a turbulent boundary layer in 2d'},
#+ {'url': 'chkpt.0022-z_face-vz-gray.png',
#+  'src': 'chkpt.0022-z_face-vz-gray.png'}]}

---

Mara is a tubulence code for astrophysical hydrodynamics and relativistic MHD.

It is the result of its author's graduate study at New York Univerity between
Summer 2008 and Spring 2013, under Andrew MacFadyen.

_Overview of capabilities_:

+ robust second-order Godunov solvers for relativistic magnetohydrodynamics
+ fifth order conservative finite-difference `WENO` schemes for the relativistic
  hydrodynamic equations
+ well-tested modules for massively parallel execution on world-class computing
  clusters, such as [Pleiades][1], efficient use of the parallel HDF5 routines
+ turbulence driving modules
+ support for arbitrary, tabulated equations of state
+ static and self-gravity coupled to the Euler equations
+ static mesh refinement, currently in 1d only

_In progress_:

+ static and adaptive mesh refinement for 2d and 3d
+ momentum and energy conserving self-gravity in 3d


[1]: http://www.nas.nasa.gov/hecc/resources/pleiades.html

### Contributors
#+{'figures':
#+[{'src': 'SRG-2048.0017-x_face-log10vort_v-afmhot.png',
#+  'url': 'SRG-2048.0017-x_face-log10vort_v-afmhot.png'}]}

---

Mara is maintained by

<address>
	<strong>Jonathan Zrake</strong><br>
    <a href="mailto:#">zrake .at. nyu .dot. edu</a><br>
    New York University<br>
    4 Washington Place<br>
    New York, NY 10003<br>
</address>

Valuable code and intellectual contributions have been made by

+ Andrew MacFadyen
+ Weiqun Zhang
+ Paul Duffell
+ Bez Laderman
+ Micha Gorelick

I also owe a great deal to

+ conversations with James Stone, and to reading parts of the [Athena code][1]
+ talks with [Tom Abel][2] at Stanford, in particular concerning mesh
  refinement
+ the [Pluto code][3], and its author Andrea Mignone
+ the excellent site [cococubed][4] maintained by Frank Timmes, in particular
  for implementing the equation of state
+ materials made available by [Christian Ott][5]
+ folks at the [Pleiades][6] help desk, who have patiently helped me to tune the
  code's parallelism and disk interaction

[1]: https://trac.princeton.edu/Athena
[2]: https://physics.stanford.edu/people/faculty/tom-abel
[3]: http://plutocode.ph.unito.it
[4]: http://cococubed.asu.edu
[5]: http://www.tapir.caltech.edu/~cott
[6]: http://www.nas.nasa.gov/hecc/resources/pleiades.html

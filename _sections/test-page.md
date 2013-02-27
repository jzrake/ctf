#+ dict(title='Mara: Test Problems', tagline='test problems')
#+ # subnav=SubNavigation([('test-1d', '1d Problems'),
#+ #                       ('test-2d', '2d Problems')], label='Contents'))

## One dimensional smooth problems

### Density wave

#### Problem description
#+ SingleImage(src='DensityWave.png',
#+             url='DensityWave.pdf')

This problem sets up a smooth density wave advecting with a uniform velocity
and periodic boundaries. This is probably the easiest problem for a code
which solves the Euler equations, as all the non-linear features are absent,
and the behavior is only that of the scalar conservation law $\dot \rho + u_0
\rho' = 0$. The problem is described by the following parameters:

+ sound speed: $c_s$
+ advection velocity: $u_0$
+ background density: $\rho_0$
+ wave amplitude: $\rho_1$
+ background pressure: $p_0 = \rho_0 c_s^2 / \Gamma$
+ wave-number: $k_0 = 2\pi / L$

$$ \left(\begin{array}{c} \rho(x,t) \\ p(x,t) \\ u(x,t) \end{array}\right) =
\left(\begin{array}{c} \rho_0 + \rho_1 \cos{k_0 (x - u_0 t)} \\ p_0 \\ 0 \\
\end{array}\right) $$

Because of its smoothness and simplicity, this problem is ideal for gauging
the convergence order of a scheme. Any scheme claiming to converge at order
$n$ must absolutely pass this test before any other. For the convergence
test, the problem is run for half of the domain-crossing time at resolutions
$2^n$ for $n$ between 3 and 10. Five different reference schemes are used:

1. `HLLC-PLM-MUSCL`
2. `HLLC-PLM-RK3`
3. `HLLC-WENO5-RK3`
4. `CHAR-WENO5-RK3`
5. `CHAR-WENO5-RK4`


#### Performance of the schemes
#+SingleImage(src='DensityWave-L1.png',
#+            url='DensityWave-L1.pdf',
#+        caption='L1 Error of various schemes')

For the formally 5th order `WENO5` reconstruction schemes, it is evident that
the truncation error introduced by the 3rd order temporal update comes to
dominate the error by resolution $64^3$. At higher resolutions, the classic 4th
order Runge-Kutta update preserves 5th order convergence. Since this problem
does not involve any non-linear waves, it is not surprising that the Godunov
scheme `HLLC-WENO5-RK3` accomplishes the same convergence order as the
characteristic decomposition scheme `CHAR-WENO5-RK3`.


| $N$  |HLLC-PLM-MUSCL|HLLC-PLM-RK3|HLLC-WENO5-RK3|CHAR-WENO5-RK3|CHAR-WENO5-RK4 |
|------|--------------|------------|--------------|--------------|---------------|
|   16 |2.07 |2.27 |4.52 |4.46 |4.53 |
|   32 |2.23 |2.04 |4.42 |4.69 |4.95 |
|   64 |2.25 |2.08 |3.97 |4.33 |5.04 |
|  128 |2.27 |2.00 |3.40 |3.66 |5.04 |
|  256 |2.36 |2.01 |3.13 |3.23 |5.13 |
|  512 |2.28 |1.98 |3.04 |3.06 |5.28 |
| 1024 |2.35 |2.00 |3.01 |3.02 |5.44 |


---

### Sound Wave


#### Problem description

This problem sets up a smooth and small-amplitude perturbation to create a
left-going sound wave. It is only slightly more challenging as a convergence
test than the density wave. The problem is described by the following
parameters:

+ sound speed: $c_s$
+ background density: $\rho_0$
+ background pressure: $p_0 = \rho_0 c_s^2 / \Gamma$
+ density wave amplitude: $\rho_1 = 10^{-6} \rho_0$
+ pressure wave amplitude: $p_1 = c_s^2 \rho_1$
+ velocity wave amplitude: $u_1 = c_s \rho_1 / \rho_0$
+ wave-number: $k_0 = 8\pi / L$
+ frequence: $\omega_0 = c_s k_0$

$$ \left(\begin{array}{c} \rho(x,t) \\ p(x,t) \\ u(x,t) \end{array}\right) =
\left(\begin{array}{c} \rho_0 + \rho_1 \cos{(k_0 x - \omega_0 t)} \\ p_0 + p_1
\cos{(k_0 x - \omega_0 t)} \\ u_0 + u_1 \cos{(k_0 x - \omega_0 t)} \\
\end{array}\right) $$

#### Performance of the schemes
#+SingleImage(src='SoundWave-L1.png',
#+            url='SoundWave-L1.pdf',
#+        caption='L1 Error of various schemes')

The need for higher temporal order is highlighted by this problem. While RK4
accomplishes 5th order convergence, there is no resolution for which RK3
converges faster than 3rd order. Notice that the error saturates at $10^{-12}$,
which is slightly larger than the roundoff error. This is because the solution
is only exact to linear order, so errors introduced at second order in $\rho_1$
are at around the $10^{-12}$ level.

| $N$  |HLLC-PLM-MUSCL|HLLC-PLM-RK3|HLLC-WENO5-RK3|CHAR-WENO5-RK3|CHAR-WENO5-RK4 |
|------|--------------|------------|--------------|--------------|---------------|
|   16 |  0.81 |0.36 |  0.61 |  0.61 |  0.78|
|   32 |  2.14 |0.47 |  2.40 |  2.40 |  3.45 |
|   64 |  2.06 |2.29 |  3.39 |  3.39 |  4.85 |
|  128 |  1.75 |1.72 |  3.18 |  3.18 |  4.71 |
|  256 |  1.95 |1.69 |  3.07 |  3.07 |  3.71 |
|  512 |  2.09 |1.85 |  3.12 |  3.12 |  0.23 |
| 1024 |  2.13 |1.95 |  2.08 |  2.08 |  0.01 |

---

## One dimensional two-state problems

The two-state problems in this section demonstrate the `Mara`'s shock-capturing
capability. Each problem sets up a piecewise-constant initial data, defined by a
different pressure, density, and velocity to either side of the midpoint of the
domain at $x=0.5$. Around discontinuities, the code will not converge at higher
order. Instead, we expect better schemes to introduce a smaller degree of
numerical dissipation around the discontinuity.

These problems show the solution as generated by a reference scheme,
`HLLC-PLM-RK3` on 128 grid points, compared against the exact
solution. [The exact Riemann solver for the Euler equation][1] is implemented
according to the algorithm of Toro (1997).

[1]: https://github.com/jzrake/ctf/blob/fd5d576ad7af8cd3d5f0e295d274e7390f039ff6/fish/riemann.c#L352


### Shocktube 1
#+ SingleImage(src='Shocktube1.png', url='Shocktube1.pdf')

This is the classic Brio-Wu shocktube problem. The solution consists of a
right-going shock wave and contact discontinuity, and a left-going rarefaction
wave. This problem is not challenging, but is excellent for debugging new
schemes because the correct solution is so familiar. Things to look for:

+ The presence of post-shock oscillations- there should be none
+ Smearing of the contact- second order (and higher) schemes will capture the
  contact in 4-6 zones

|        |  $x<0.5$ |  $x>0.5$ |
|--------|----------|----------|
| $\rho$ | 1.000 | 0.125 |
| $p$    | 1.000 | 0.100 |
| $v_x$  | 0.000 | 0.000 |
| $v_y$  | 0.000 | 0.000 |
| $v_z$  | 0.000 | 0.000 |

---

### Shocktube 2
#+ SingleImage(src='Shocktube2.png', url='Shocktube2.pdf')

Double rarefeaction wave: in this setup, the density and pressure are the same
across the midline. The velocity is discontinuous, with both halves of the
domainmoving away from midline. This problem is helpful for finding asymmetries
in the code. It should be run on both an even and odd number of zones, to
confirm that the solution is always left-right symmetric.

|        |   $x<0.5$ |  $x>0.5$ |
|--------|-----------|----------|
| $\rho$ |  1.000 | 1.000 |
   | $p$ |  0.400 | 0.400 |
| $v_x$  | -2.000 | 2.000 |
| $v_y$  |  0.000 | 0.000 |
| $v_z$  |  0.000 | 0.000 |

---

### Shocktube 3
#+ SingleImage(src='Shocktube3.png', url='Shocktube3.pdf')

This setup creates a strong right-going shock wave by initializing a pressure
jump of 5 orders of magnitude. The quality of the solution can be gauged by the
percentage of the true maximal density behind the shock. In the example solution
given, with $128$ zones, there are two zones which reach roughly 80% the
maximal density.

|        |     $x<0.5$ |  $x>0.5$ |
|--------|-------------|----------|
| $\rho$ |    1.000 | 1.000 |
| $p$    |     1000 | 0.010 |
| $v_x$  |    0.000 | 0.000 |
| $v_y$  |    0.000 | 0.000 |
| $v_z$  |    0.000 | 0.000 |

---

### Shocktube 4
#+ SingleImage(src='Shocktube4.png', url='Shocktube4.pdf')

This setup creates a strong left-going shock wave by initializing a pressure
jump of 4 orders of magnitude.

|        |  $x<0.5$ |    $x>0.5$ |
|--------|----------|------------|
| $\rho$ | 1.000 |   1.000 |
| $p$    | 0.010 |   100.0 |
| $v_x$  | 0.000 |   0.000 |
| $v_y$  | 0.000 |   0.000 |
| $v_z$  | 0.000 |   0.000 |

---

### Shocktube 5
#+ SingleImage(src='Shocktube5.png', url='Shocktube5.pdf')

This setup yields a two-shock solution by creating a colliding supersonic
flow. The velocities are adjusted so that the left shock is stationary.

|        |    $x<0.5$ |   $x>0.5$ |
|--------|------------|-----------|
| $\rho$ |   5.999240 |  5.999240 |
| $p$    | 460.89400  | 46.0950   |
| $v_x$  |  19.597500 | -6.196330 |
| $v_y$  |   0.000    |  0.000    |
| $v_z$  |   0.000    |  0.000    |

---

### Isolated Contact Wave
#+ SingleImage(src='ContactWave.png', url='ContactWave.pdf')

This problem consists of only a stationary jump in the density. Since the Euler
equations do not model any diffusive processes, the jump should remain
undisturbed. However, the `HLL` riemann solver will, in fact create anomalous
diffusion across the contact. The contact wave problem is useful at verifying an
`HLLC` solver, which uses a two-state rather than a one-state approximation to
the intermediate state in order to capture the contact discontinuity.

|        |  $x<0.5$ |  $x>0.5$ |
|--------|----------|----------|
| $\rho$ | 1.000 | 0.100 |
| $p$    | 1.000 | 1.000 |
| $v_x$  | 0.000 | 0.000 |
| $v_y$  | 0.700 | 0.700 |
| $v_z$  | 0.200 | 0.200 |

---


## Two-dimensional test problems

### Kelvin-Helmholtz instability
#+ SingleImage(src='KH-hllc-plm-rk3.png',
#+             url='KH-hllc-plm-rk3.png',
#+         caption='Performance at various resolutions for the '
#+                 '<code>HLLC-PLM-RK3</code> scheme')

This problem uses a smooth shearing profile to capture the linear growth rate of
Kelvin-Helmholtz instability. The vertical velocity is given a sinusoidal
perturbation with 4 wave-lengths over the domain. The resulting flow, when
properly resolved, contains a single vortex for each wave-length of the
perturbation. The domain is $[0,L]^2$.

+ background pressure: $p_0 = 2.5$
+ outer density: $\rho_1 = 1.0$
+ inner density: $\rho_2 = 2.0$
+ outer velocity: $u_1 = -0.5$
+ inner velocity: $u_2 =  0.5$
+ shearing layer width: $\delta = 0.035$
+ perturbation amplitude: $w_0 = 10^{-2}$

 $$ \left(\begin{array}{c} \rho \\ p \\ u \\ v \\ \end{array}\right) =
 \left(\begin{array}{c} \frac{\rho_2 - \rho_1}{2} (\tanh{\frac{y-L/4}{\delta}} -
 \tanh{\frac{y-3L/4}{\delta}}) + \rho_1 \\ p_0 \\ \frac{u_2 - u_1}{2}
 (\tanh{\frac{y-L/4}{\delta}} - \tanh{\frac{y-3L/4}{\delta}} - 1) \\ w_0 \sin(4
 \pi x) \\ \end{array}\right) $$

---

### Implosion with reflecting walls

#### Problem description
#+ Movie('Implosion2d-HLLC-PLM-MUSCL.mp4')

This setup involves similar conditions as in the `Shocktube1` problem, except
that the discontinuity is placed at a $45^{\circ}$ angle on a 2d grid, and the
boundary conditions are reflecting. The higher pressure region (above the
diagonal) causes a shock to propagate toward the lower left corner, which then
reflects off the other walls multiple times. The interaction of the shock with
various contact discontinuities helps diagnose the ability of various schemes to
limit artificial dissipation of contacts. This test is also a stringent test of
the code's $x-y$ symmetry. If the problem is run sufficiently long, any bug
which breaks the reflectional symmetry will show up by skewing the features
above or below the diagonal.

$$ \left(\begin{array}{c} \rho \\ p \\ u \\ v\end{array}\right) =
\begin{cases}
\left(\begin{array}{c} 1.000 \\ 1.000 \\ 0 \\ 0 \end{array}\right) & x + y > L/2 \\
\left(\begin{array}{c} 0.125 \\ 0.140 \\ 0 \\ 0 \end{array}\right) & \text{otherwise}
\end{cases} $$

This problem is also documented on the [Athena test page][1].

[1]: http://www.astro.princeton.edu/~jstone/Athena/tests/implode/Implode.html

#### Performance of the schemes
#+ TabbedImages([dict(src='HLLC-PLM-MUSCL.0480.png', alt='HLLC-PLM-MUSCL'),
#+               dict(src='HLLC-PLM-RK3.0480.png', alt='HLLC-PLM-RK3'),
#+               dict(src='CHAR-WENO5-RK3.0480.png', alt='CHAR-WENO5-RK3'),
#+               dict(src='CHAR-WENO5SZ-RK3.0480.png', alt='CHAR-WENO5SZ-RK3')])

The suite of schemes chosen for this problem is: two which use the HLLC
approximate Riemann solver with a 2nd order Godunov scheme based on the
piecewise linear method, and two characteristic-wise WENO schemes.

Of the Godunov schemes, one is `HLLC-PLM-MUSCL`, which is 2nd order in time and
uses an _unsplit_ integration scheme with corner-transport while the other,
`HLLC-PLM-RK3`, uses an unplit method-of-lines approach with a three-stage,
third order single-register (low storage) Runge-Kutta
integration. `HLLC-PLM-MUSCL` is particuarly robust for relatistic MHD
turbulence problems with moderate $\sigma$. One reason for this is apparent from
the figures at the right: the `MUSCL` scheme is more diffusive. Still, as shown
in the smooth problems, that scheme is formally 2nd order accurate.

The other two schemes are dimensionally-split, characteristic-wise conservative
finite-difference schemes based on the method of lines and the same 3rd order
Runge-Kutta integrator. These schemes are both formally 5th order accurate
(although with RK3, they transition to 3rd order at around 128 zones) but differ
in the _smoothness indicators_ they utilize. `CHAR-WENO5-RK3` uses the older
smoothness indicators ($IS_k$) of Jiang and Shu (1996), whereas
`CHAR-WENO5SZ-RK3` uses the improved $IS_k$ of Zhang and Shen (2010). The
performance of the updated $IS_k$ is so superior to the older ones that these
schemes are the most, and least diffusive of those tested here.

---

## Self-gravity

### Sound wave with gravity

#### Problem description
#+ SingleImage(src='SoundWaveGrav.png',
#+             url='SoundWaveGrav.pdf',
#+         caption='Gravitating sound wave evolved for slightly less than '
#+                 'one sound crossing time')

The presence of gravity adds a correction at linear order to the acoustic
dispersion relation. The linearized Euler equations for the mass and momentum
conservation, respectively read

$$ \rho_0 \dot u + p' + \rho_0 \phi' = 0 \\ \dot \rho + \rho_0 u' = 0 $$

Together with the linearized equation of state $p' = c_s^2 \rho'$, where $c_s^2
= \Gamma p_0 / \rho_0$ is the sound speed, and the Poisson equatio, $\phi'' =
\rho - \rho_0$ this system can be written as

$$ c_s^2 \rho'' - \ddot \rho + \rho_0 (\rho - \rho_0) = 0 $$

Assuming $\rho -\rho_0 \propto \cos{kx - \omega t}$, the dispersion relation
falls out:

$$ \omega^2 = c_s^2 k^2 - \rho_0 $$

Here, it may be helpful to re-dimensionalize by replacing $4 \pi G$, which until
now has been set to unity. Then the dispersion relation in terms of the Jeans
$\tau_J = 1/\sqrt{G \rho}$ time is $ \omega^2 = c_s^2 k^2 - 4\pi / \tau^2 $. The
wave is marginally stable against gravitational collapse when the right hand
side is zero. This occurs for density perturbations at the Jeans length
$\lambda_J = c_s \sqrt{\pi / G \rho_0}$.

	
### Collapse of a cold dust cloud

#### Problem description
#+ SingleImage(src='Collapse1d-rhox.png', url='Collapse1d-rhox.pdf')

This problem tests the collapse of a cold dust cloud in one dimension. The
initial density field is a uniform-density region of width $\delta$ with a
low-density "atmosphere":

$$ \rho(x) = \begin{cases} \rho_0 & |x|<\delta/2 \\ \rho_\rm{atm} &
\text{otherwise} \end{cases} $$

The problem domain is $[-L/2, L/2]$ with periodic boundaries so that the total
mass $M = \delta \rho_0$ is a constant. The exact potential can easily be
written down for this density field. In periodic boundaries, solutions only
exist when the net charge is zero, so we solve the Poisson equation

$$ \nabla^2 \phi = \rho - \bar{\rho} $$

for deviations away from the mean density $\bar{\rho} = M/L$. Letting
$\phi(\pm \delta/2) = 0$, we can immediately write down the potential.

$$ \phi(x) = \begin{cases} -\frac{\rho_0}{2} (x+\delta/2)(x+3\delta/2+L) & x <
-\delta/2 \\ -\frac{\rho_0}{2} (x-\delta/2)(x-3\delta/2-L) & x > \delta/2 \\
\frac{\rho_0(1-\delta/L)}{2} (x-\delta/2)(x+\delta/2) & \text{otherwise} \\
\end{cases} $$

Under the influence of gravity, all fluid elements experience an acceleration
$a(x) = -\phi'(x) = -\rho_0 (1 - \delta/L) x$. As the dust cloud collapses, the
density $\rho_0(t)$ increases while $\delta(t)$ decreases so as to hold $M$
fixed. By following the location of the density jump $\delta(t)/2$, we can write
down the following ODE for $\delta(t)$.

$$ \ddot{\delta} + M(1 - \frac{\delta}{L}) = 0 $$

which has the solution $\delta(t) = L - (L - \delta_0) \cosh(t/\tau)$, where
$\tau = (M/L)^{-1/2}$ is the Jeans time.

#### Code performance
#+ SingleImage(src='Collapse1d-rhot.png', url='Collapse1d-rhot.pdf')

The code performs well on this problem, provided the atmosphere $\rho_\rm{atm}$
is sufficiently large, $10^{-6} \times \rho_0$ works fine. The density profile
remains reasonably square, but with some slight artifact at the density jump. In
order to compare with the exact solution, the maximum density on the grid is
recorded at each time step, and compared against the true solution

$$ \rho_0(t) = \frac{M/L}{1 - (1-\delta_0/L) \cosh(t/\tau)} $$

---

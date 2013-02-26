#+
#+{'figures': [{'src': 'DensityWave-L1.png',
#+              'url': 'DensityWave-L1.pdf',
#+              'caption': 'L1 Error'}]}
#+

### Problem description

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


### Performance of the schemes

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


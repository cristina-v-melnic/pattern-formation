# pattern-formation

## About

Project that models the spontaneous formation of patterns in reaction-diffusion system 
based on Alan Turing's paper of 1952 "The Chemical Basis of Morphogenesis" [1], where,
the development of distinct biological features is modelled as a 
coupled reaction-diffusion process between two morphogens $u$ and $v$, i.e.,

$$
\partial_{t} u = \delta_{1} \Delta u + f(u,v),
$$

$$
\partial_{t} v = \delta_{2} \Delta v + g(u,v).
$$

Here, we solved numerically the above system of differential equations using the
Euler and the Crank-Nicolson integration rules and found a better performance
in the latter case, as expected due to its implicit nature and a higher convergence order.
The computational expense of an implicit integrator (solving a nonlinear fixed-point problem)
is, in this case, compensated with by the property of unconditional stabillity.
This allows integration with large time steps $\Delta t > 1$, which considerably reduces the computational time. 
The computational properties were studied using models from literature for abstract pattern generation [2], as
well as formation of animal fur patterns [3-4], which were specified by the values of parameters $\delta_{1}$ and $\delta_{2}$
and functional dependencies $f(u,v)$ and $g(u,v)$. 

## Credit
This work was done as a part of a semester project in MSc Computational Sciences at the Freie Universitaet Berlin in 2019. The project was initiated and supervised by Dr Gottfried Hastermann and performed by a group of 3 people: myself, Ece Sanin and Kevin Cyriac Edampurath. 


## References

[1]. A. M. Turing, “The chemical basis of morphogenesis", Philosophical Transactions of the Royal Society of London. Series B, Biological Sciences, vol. 237, pp. 37{72, Aug. 1952.} 

[2]. J. E. Pearson, “Complex Patterns in a Simple System,” Science, vol. 261, pp. 189–192, July 1993. 		

[3]. A. J. Koch and H. Meinhardt “Biological Pattern Formation : from Basic Mechanisms to Complex Structures”, Rev. Modern Physics 66, 1481-1507 (1994)

[4]. R. T. Liu, S. S. Liaw, and P. K. Maini, "Two-stage Turing model for generating pigment patterns on the leopard and the jaguar," Physical Review E, vol. 74, p. 011914, July 2006.

## Contents

- *solvers.py* - Convection-diffusion solver classes that implement numerical integration of coupled convection-diffusion
equations to get the solutions $u(t,x,y)$ and $v(t,x,y)$;
- *models.py* - Classes that specify the type of parameters and functional dependencies between the 2 unknowns u and v;
- *performance.py* - Benchmarking function.


### TO DO:
- [ ] Add description to classes and functions.
- [ ] Update the test.

# Manual for the augmented MERCURY integrator


This version of the [OrbFit](http://adams.dm.unipi.it/orbfit/) package contained a modified version of the orbit9 integrator. The new features included in this package are specifically designed for the long-term propagation of small solar system objects, i.e. asteroids.

More specifically, the Yarkovsky effect is added to the gravitational vector field, and the spin-axis evolution due to the YORP effect is integrated together with the orbital dynamics. Technical details of the equations and of the implementation can be found in the paper

"Mercury and OrbFit packages for numerical integration of planetary systems: implementation of the Yarkovsky and YORP effects", M. Fenucci and B. Novaković

If you publish results using this version of the integrator, please reference the package using the above paper. 

In this readme you can find instructions on how to use the additional routines provided in this version of the integrator.


## Authors 
- [Marco Fenucci](http://adams.dm.unipi.it/~fenucci/index.html), Department of Astronomy, Faculty of Mathematics, University of Belgrade (<marco_fenucci@matf.bg.ac.rs>) 
- [Bojan Novaković](http://poincare.matf.bg.ac.rs/~bojan/index_e.html), Department of Astronomy, Faculty of Mathematics, University of Belgrade (<bojan@matf.bg.ac.rs>) 

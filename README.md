# Manual for the augmented MERCURY integrator


This software is an augmented version of the MERCURY integrator, that was  originally developed by John E. Chambers. The new features included in this package are specifically designed for the propagation of small solar system objects, i.e. asteroids.

More specifically, the Yarkovsky effect is added to the gravitational vector field, and the spin-axis evolution due to the YORP effect is integrated together with the orbital dynamics. Technical details of the equations and of the implementation can be found in the paper

" ", M. Fenucci and B. Novaković, Serbian Astronomical Journal

If you publish results using this version of the integrator, please reference the package using the above paper. 

In this readme you can find instructions on how to use the additional routines provided in this version of the integrator. To prepare the initial conditions for planets and small objects, please refer to the original manual written by John E. Chambers. 


## Authors 
- [Marco Fenucci](http://adams.dm.unipi.it/~fenucci/index.html), Department of Astronomy, Faculty of Mathematics, University of Belgrade (<marco_fenucci@matf.bg.ac.rs>) 
- [Bojan Novaković](http://poincare.matf.bg.ac.rs/~bojan/index_e.html), Department of Astronomy, Faculty of Mathematics, University of Belgrade (<bojan@matf.bg.ac.rs>) 

# Manual for the augmented orbit9 integrator

The [OrbFit](http://adams.dm.unipi.it/orbfit/) software was originally developed by Prof. Andrea Milani Comparetti (University of Pisa), and it is currently maintained and developed by the OrbFit consortium. OrbFit offers many tools to the user, including programs for orbit determination, orbit propagation, and ephemeris computation.

This version of the [OrbFit](http://adams.dm.unipi.it/orbfit/) package contains a modified version of the orbit9 integrator, a program that permits to propagate the dynamics of asteroids. The new features included in this package are specifically designed for the long-term propagation of small solar system objects, i.e. asteroids.

More specifically, the Yarkovsky effect is added to the gravitational vector field, and the spin-axis evolution due to the YORP effect is integrated together with the orbital dynamics. Technical details of the equations and of the implementation can be found in the paper

"Mercury and OrbFit packages for numerical integration of planetary systems: implementation of the Yarkovsky and YORP effects", M. Fenucci and B. Novaković

If you publish results using this version of the integrator, please refer to the package using the above paper. 

In this readme you can find instructions on how to use the additional routines provided in this version of the integrator.


## Authors 
- [Marco Fenucci](http://adams.dm.unipi.it/~fenucci/index.html), Department of Astronomy, Faculty of Mathematics, University of Belgrade (<marco_fenucci@matf.bg.ac.rs>) 
- [Bojan Novaković](http://poincare.matf.bg.ac.rs/~bojan/index_e.html), Department of Astronomy, Faculty of Mathematics, University of Belgrade (<bojan@matf.bg.ac.rs>) 

## Compilation

The compilation is the same as for the OrbFit version sdistributed by the OrbFit Consortium. The distribution comes with a bash script called config, that permits to choose the compilation settings for several FORTRAN compilers. By running the script without any option, the user will receive a help message. The most popular FORTRAN compilers are the GNU gfortran, and the INTEL ifort, that can be chosen with the script. 

We suggest the final user to choose the optimization flags, by

        ./config -O gfortran
        
At this point, the code can be compiled by typing 

        make

**Note.** This operation may take several minutes to complete.

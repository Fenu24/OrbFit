2014/08/07
How to create new JPL ephemerides.

1) Copy from the web the useful file (something similar to ascp0*).
2) Pay attention to the README files given by the JPL. asc2eph.f90 and
testeph.f90 are updated (last modified: 2013, 09). If you need new
files you have to download them from the web site. Usually JPL use f77
files, which means that you have to convert the files. Use the routine
convert.x in util90.
3) make
4) Change the Makefile in order to have what you need, e.g. 
   input: 
        cat header.431_572 ascp00000.431 ascp01000.431 ascp02000.431 ascp03000.431  > input.431 
 
   is for DE431 from 0 to 4000. If you want you can change the time interval.
5) make ephemerides

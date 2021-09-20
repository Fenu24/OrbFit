For your convenience we have implemented a graphical interface
"OrbFitSoft" in Orbfit (Unix/Linux only), which is intended to
facilitate usage of the entire system. You can run all the programs
contained in the package and make use of the complete, online help,
available even while other programs are running.

We assume that you have a Tcl/Tk package installed, and that 'wish'
shell is in '/usr/bin/'. If this is not the case you can either set an
appropriate link, or you can edit all the executables in the 'doc/'
subdirectory of Orbfit (OrbFitSoft, HELP, ORBFIT, FITOBS, BINEPH) and
change the first line in each of them to reflect the location of
'wish' on your system.

This application is developed using version 8.0 of 'wish'. To learn
more on Tcl/Tk consult either an online 'tclhelp' facility, available
in '/usr/bin/', or the book "Tcl and the Tk Toolkit" by J.K.Ousterhout
(Adison-Wesley Publ.Comp.).

To use "OrbFitSoft" just invoke it from the command line in your
working directory. You should, however, first set a couple of
links in your working directory (links are already provided in all 
the 'tests/<name>' subdirectories); if, for example, the working
directory is 'home/zoran/workdir', and Orbfit is installed in
'/home/zoran/orbfit', then you should:

% cd /home/zoran/workdir
% ln -s ../orbfit/doc/OrbFitSoft
% ln -s ../orbfit/src/include/doclib.h

Your graphical interface is now fully operational and you can type

% OrbFitSoft

to invoke it. If you lose contents of the windows when you drag them
from one place on the screen to another or overlap one window with
another then you should not run OrbFitSoft in the background.

Usage of the interface is simple and essentially self-explanatory (for
additional information see also "Help on Help" in menu "ABOUT" of the
help facility).

Note, however, that the parent 'xterm' window from which you invoke
OrbFitSoft will be grabbed by the application and not released until
you exit from the graphical interface. Thus, all the screen output of
the current program is going to appear in this window.  Also you will
not get the prompt back in this window after launching any of the
jobs, even when the current job ends (thus do not wait for the prompt
to learn that the job is over, but there will be an execution time
report), until you release it by exiting the interface.  Still, you
can relaunch the current job, or launch any other job by simply
clicking on the corresponding menu item in the top-level "AstOrb"
window. If you need to do something else while the application is
running, you should open another 'xterm' window.

Help is available even while other programs are running. Other
programs (orbfit, fitobs, bineph) should preferably run one at a time,
as their outputs would otherwise be all mixed up in the same window.

The automatic placement of the interface windows is optimized for a
17-inch screen, with 1024 x 768 resolution. If you have problems with
the positioning of the window and/or the menus (submenus getting out
of the screen), just edit the last line of the corresponding file (see
above for the list of executables), which shoud be something like:

wm geometry . +160+20

and change the numbers (in pixels on the display counting from the 
upper-left corner) until you find an appropriate position for the window 
on your screen.
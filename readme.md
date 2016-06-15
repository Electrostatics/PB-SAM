#This repository has been deprecated.
##Building PB-SAM

This is the directory containing source code and examples for the semi-analytical solution of
the linearized Poisson Boltzmann equation, as described in Yap and Head-Gordon 2010/2012.

The src directory contains source code.  To make, open the makefile, change any
libs or includes as needed and make as follows:

	For make sphere, cd makesphere, then
		make makesphere
	For the rest, cd ../, then
		make findContacts
		make pbsam_30
	## Before making any bd, edit the file expansion.h: #define N_POLES 10 
		make bdnam_10
		make bdm_10

The bin directory stores the executable, called makesphere, findContacts, pbsam, bdnam, bd.

The test directory contains a series of submit scripts and required inputs for each
step of the PB-SAM code.  Each is described in the README in the test directory.

The doc directory contains doxygen created documentation for the source code.  It was
created by doxygen version 1.8.1.  To load as an html, open the annotated.html file in
a browser.

## Building PB-SAM on OSX

* Install gcc (homebrew will do this)
* Edit the src/makefile:
 * +CC = g++-5
 * +LINK =  -L$(MKL_HOME)/lib -lmkl_intel_lp64 -lmkl_core -ldl -lpthread -lm

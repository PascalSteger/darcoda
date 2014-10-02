GravLite: Mass modelling tool for spherical and disk-like structures
====================================================================


Introduction
------------

GravLite is a tool to determine the mass distribution in disc-like and spherical
systems. It takes as input a tracer density distribution, a line-of-sight
velocity dispersion, and possibly the velocity's fourth moment as a function of
radius. It then generates a highly-dimensional parameter space for tracer
density, overall density distribution, and a velocity anisotropy profile in bins
and samples it with MultiNest, a Monte Carlo method. The 1-dimensional Jeans
equations for these systems is then used to calculate a goodness-of-fit from the
surface density and velocity dispersion along the line of sight. The accepted
models can be visualized in a later step.


Installation
------------

Following packets need to be installed on your system:
 * python3, ipython3
 * matplotlib/pylab
 * scipy
 * ipdb

Unzip the file gravlite.zip, and run

> ipython3 gravlite.py;


Parameter files: Main configuration
-----------------------------------

Sample parameter files are stored in the subfolders disc/ and sphere/. The file
./gl_params.py is a soft link to one of them. Following mass modelling methods
have been implemented so far:

 * walk: spherical Walker mock data from the Gaia challenge catalogue, 2 populations
 * gaia: spherical mock data from the Gaia challenge catalogue, 1 population
 * hern: spherical mock data taken from a Hernquist profile
 * discmock: disk-like mock data, generated on the fly
 * discsim: disk-like mock data, from a simulation by S. Garbari

Further Documentation
---------------------

Additional documentation on profile representation, priors, runtime parameters
and more can be found in the paper draft folder, and in the code itself.

A HTML documentation of all files and functions can be generated using

> doxygen Doxyfile

and then browsing to doc/html/.

Bugs
----

The code is under constant development, and might show bugs. Feel free to branch
the code, correct them, and send a patch!


October 2014,
Pascal Steger
psteger@phys.ethz.ch
http://www.n.ethz.ch/~psteger

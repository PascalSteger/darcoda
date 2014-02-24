Gravlite: Mass modelling tool for spherical and disk-like structures
====================================================================


Introduction
------------

Gravlite is a tool to determine the mass distribution in disc-like and spherical systems. It takes as input a tracer density distribution, a line-of-sight velocity dispersion, and possibly the velocity's fourth moment as a function of radius. It then generates a highly-dimensional parameter space for tracer density, overall density distribution, and a velocity anisotropy profile ands samples it with MultiNest, a Monte Carlo method.


Installation
------------

Following packets need to be installed on your system:
 * python3
 * matplotlib/pylab
 * scipy

Unzip the file gravlite.zip, and run

python3 gravlite.py;


Parameter files: Main configuration
-----------------------------------

Sample parameter files for several possible scenarios are stored in the subfolder ./params. The file gl_params.py is a soft link to one of them (gl_params_gaia.py is constantly kept updated, the others might need modifications). Following mass modelling methods have been implemented so far:

 * gl_params_walker.py: spherical Walker mock data from the Gaia challenge catalogue, 2 populations
 * gl_params_gaia.py: spherical mock data from the Gaia challenge catalogue, 1 population
 * gl_params_hern.py: spherical mock data taken from a Hernquist profile
 * gl_params_simple.py: disk-like mock data, generated on the fly
 * gl_params_sim.py: disk-like mock data, from a simulation by S. Garbari

Further Documentation
---------------------

Additional documentation on profile representation, priors, runtime parameters and more can be found in the paper draft folder, and in the code itself. Start in gravlite.py, and have fun!

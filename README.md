GravImage: Mass modelling tool for spherical and disk-like structures
=========]===========================================================

Introduction
------------

GravImage is a tool to determine the mass distribution in
one-dimensional disc-like or spherical systems. It takes as input a
tracer density distribution, a line-of-sight velocity dispersion, and
possibly the velocity's fourth moment as a function of radius. It then
generates a highly-dimensional parameter space for tracer density,
overall density distribution, and a velocity anisotropy profile in
bins and samples it with MultiNest, which is a specialized
high-dimensional parameter space sampling algorithm. The Jeans
equations for the systems under study is then used to calculate a
goodness-of-fit from the surface density and velocity dispersion along
the line of sight. The accepted models are visualized in a later
step.


Installation
------------

Following packages need to be installed on your system:
 * python3
 * matplotlib/pylab
 * scipy
 * ipdb, pdb

Then execute

> git clone https://github.com/PascalSteger/darcoda $DARCODA_DIR

> cd $DARCODA_DIR

or unzip the file darcoda.zip to $DARCODA_DIR . Then set the environment variables

> export PYTHONPATH=$PYTHONPATH:$DARCODA_DIR/gravimage/programs/
> export PYTHONPATH=$PYTHONPATH:$DARCODA_DIR/gravimage/programs/plotting/

Adapt the path specifications to your needs in

gl_params.py
gl_class_files.py
import_path.py

and run

> python3 gravimage.py



Parameter files: Main configuration
-----------------------------------

Sample parameter files are stored in the subfolders disc/ and sphere/. The file
./gl_params.py is a soft link to one of them. Following mass modelling methods
have been implemented so far:

 * hern: spherical mock data taken from a Hernquist profile
 * gaia: spherical mock data from the Gaia challenge catalogue, 1 population
 * walk: spherical Walker mock data from the Gaia challenge catalogue, 2 populations
 * obs: observations of 4 dwarf spheroidals
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


December 2014,
Pascal Steger
psteger@phys.ethz.ch
http://steger.aero

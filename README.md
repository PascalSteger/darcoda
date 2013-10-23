Gravlite: Gravlite - Non-Parametric Mass Modelling Routine for Discs and Spheres

Mass modelling tool for spherical and disk-like structures

Introduction
Gravlite is a tool to determine the mass distribution in disc-like and spherical systems. It takes as input a tracer density distribution, a line-of-sight velocity dispersion, and the velocity's fourth moment as a function of radius. It then generates a highly-dimensional parameter space for tracer density, overall density distribution, and possibly a velocity anisotropy profile ands traces it with a simple Markov Chain Monte Carlo method.


Installation
Following packets need to be installed on your system:
 * python3
 * matplotlib/pylab
 * scipy

Unzip the file gravlite.zip, and run

python3 gravlite.py;


Parameter files: Main configuration

Sample parameter files for several possible scenarios are stored in the subfolder ./params. The file gl_params.py is a soft link to one of them. Following mass modelling methods have been implemented so far: 

 * gl_params_walker.py: spherical Walker mock data from the Gaia challenge catalogue, 2 populations
 * gl_params_gaia.py: spherical mock data from the Gaia challenge catalogue, 1 population
 * gl_params_hern.py: spherical mock data taken from a Hernquist profile
 * gl_params_simple.py: disk-like mock data, generated on the fly
 * gl_params_sim.py: disk-like mock data, from a simulation by S. Garbari


Init values
Initial values for all profiles and the respective stepsizes are set in gl_init.py. We use following representations:

 * overall density rho: can either be set as rho_i=rho(r_i), where r_i are the radii of bins, or the coefficients of the Legendre polynomial;
 * tracer densities nu_i = nu(r_i) in linear space or nu_i = log(nu(r_i)) in logarithmic space;
 * velocity anisotropy beta_i are incremental changes from beta(r=0) = 0.

Have fun!

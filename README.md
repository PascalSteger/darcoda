Darcoda: Codes to search for Dark Matter
========================================

Introduction
------------

GravLite is a tool to determine the mass distribution in
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


October 2014,
Pascal Steger
psteger@phys.ethz.ch
http://www.n.ethz.ch/~psteger

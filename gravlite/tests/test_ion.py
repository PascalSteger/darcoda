#!/usr/bin/python
# (c) 2013 Pascal S.P. Steger

from pylab import *

t = linspace(0.0, pi, 100)
x = cos(t)
y = sin(t)

ion()  # turn on interactive mode
figure(0)
subplot(111, autoscale_on=False, xlim=(-1.2, 1.2), ylim=(-.2, 1.2))

point = plot([x[0]], [y[0]], marker='o', mfc='r', ms=3)

for j in arange(len(t)):
    # reset x/y-data of point
    setp(point[0], data=(x[j], y[j]))
    draw() # redraw current figure

ioff() # turn off interactive mode
show()
################ 

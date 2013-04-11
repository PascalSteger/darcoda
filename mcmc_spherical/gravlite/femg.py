#!/usr/bin/python2.7
# (c) 2013 Pascal Steger, psteger@phys.ethz.ch
'''search for groups of stellar tracers in Fe/Mg space'''

import numpy
import sys
if(len(sys.argv)<2):
    print "use: femg.py [car,scl,sex,for]"
    exit(1)
    
dwarf=sys.argv[1]
dir="/home/ast/read/dark/dwarf_data/"
print dir+dwarf+"/table_merged.bin"

delim=[0,22,3,3,6,4,3,5,6,6,7,5,6,5,6,5,6]
ID=numpy.genfromtxt(dir+dwarf+"/table_merged.bin",skiprows=29,unpack=True,usecols=(0,1),delimiter=delim,dtype="string")
RAh,RAm,RAs,DEd,DEm,DEs,Vmag,VI,VHel,e_VHel,SigFe,e_SigFe,SigMg,e_SigMg,PM=numpy.genfromtxt(dir+dwarf+"/table_merged.bin",skiprows=29,unpack=True,usecols=tuple(range(2,17)),delimiter=delim,filling_values=-1)

# only use stars which are members of the dwarf
pm = (PM>0.95)*(VI<70)*(SigFe>-1)*(SigMg>-1)
print "fraction of members = ",1.0*sum(pm)/len(pm)
ID=ID[1][pm]; RAh=RAh[pm]; RAm=RAm[pm]; DEd=DEd[pm]; DEm=DEm[pm]; DEs=DEs[pm]; Vmag = Vmag[pm]; VI=VI[pm]; VHel=VHel[pm]; e_VHel=e_VHel[pm]; SigFe=SigFe[pm]; e_SigFe=e_SigFe[pm]; SigMg=SigMg[pm]; e_SigMg=e_SigMg[pm];PM=PM[pm]

from pylab import *
# from matplotlib.patches import Ellipse
# ion(); a = subplot(111, aspect='equal')

# for i in range(len(e_SigFe)):
#     e=Ellipse((SigFe[i],SigMg[i]), e_SigFe[i], e_SigMg[i], 0.)
#     e.set_clip_box(a.bbox)
#     e.set_alpha(0.01)
#     a.add_artist(e)

# xlabel(r'Fe')
# ylabel(r'Mg')
# #xlim([0,1])

import numpy as np
#import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

# the random data
x = SigFe
y = SigMg

nullfmt   = NullFormatter()         # no labels

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left+width+0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

# start with a rectangular Figure
f=figure(1, figsize=(8,8))

axScatter = plt.axes(rect_scatter)
axHistx   = plt.axes(rect_histx)
axHisty   = plt.axes(rect_histy)

# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

# the scatter plot:
#axScatter.scatter(SigFe,SigMg,'b.')
from matplotlib.patches import Ellipse
for i in range(len(e_SigFe)):
    e=Ellipse((SigFe[i],SigMg[i]), e_SigFe[i], e_SigMg[i], 0.)
    e.set_clip_box(axScatter.bbox)
    e.set_alpha(0.01)
    axScatter.add_artist(e)


# now determine nice limits by hand:
Nbins = 20
binwidth = 1./Nbins
xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
lim = ( int(xymax/binwidth) + 1) * binwidth

axScatter.set_xlim( (0, 1) )#-lim,lim
axScatter.set_ylim( (0, 1) )#-lim,lim

bins = np.arange(-lim, lim + binwidth, binwidth)
axHistx.hist(x, bins=bins)
axHisty.hist(y, bins=bins, orientation='horizontal')

axHistx.set_xlim( axScatter.get_xlim() )
axHisty.set_ylim( axScatter.get_ylim() )

savefig(dir+dwarf+"/femg.png")
ioff();show()

exit(0)

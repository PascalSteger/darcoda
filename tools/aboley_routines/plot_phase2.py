#!/usr/bin/python
'''Calculate the radial distribution from any data cube.'''

import numpy
import array
import struct
import pylab
import math
import sys
import toolKit2 as TK

# file-dependent variables.
istart=32
istop =33

RM=4
mp=1.67e-24
kb=1.3801e-16
USECENTERING=0
OUTDIR="/scratch/psteger/sim/"
DIR="./"

for i in xrange(istart,istop):
	number=repr(i).zfill(5)
	filename=DIR+"cell"+number+".dat"
	f = open(filename,'r')
	header=TK.getRamsesInfo(i,OUTDIR+"./output_")
	ncpu,ndim,levelmin,levelmax,ngridmax,nstep_coarse,boxlen,time,aexp,H0,omega_m,omega_l,omega_k,omega_b,unit_l,unit_d,unit_t=header

	norm=unit_d/mp
	polytropic_constant=( 4.*0.5**float(levelmax)*unit_l/aexp)**2/(math.pi)*6.67e-8*unit_d*(unit_t/unit_l)**2
	polytropic_constant=0.

	dd=[]
	tk=[]
	metal=[]

	l=0
	ddMax=0.
	offSet=0.
	xx=0.;yy=0.;zz=0.
	for line in f:
		values=line.split()
		dd.append(float(values[6]))
		tk.append(float(values[10]))
		tk[l]-=polytropic_constant*dd[l]**2*2./3.
		tk[l]=max(tk[l]/dd[l]*(unit_l/unit_t)**2*mp/kb,1e-1)
		metal.append(max(math.log10(float(values[11])),-9.))
		dd[l]*=norm
		if tk[l]<.1 and dd[l]<.1 and dd[l]>.05:
			xx=float(values[0])
			yy=float(values[1])
			zz=float(values[2])
		l+=1
        f.close()
	filename="zrho"+number+".png"         
	filename2="ztk"+number+".png"         
	filename3="phase"+number+".png"         

	print numpy.min(dd),numpy.max(dd)
	print numpy.min(tk),numpy.max(tk)
        print xx,yy,zz

        n=len(metal)
	r_array=[]
	g_array=[]
	b_array=[]
	zmax=max(metal)
	for i in xrange(n):
		r_array.append((metal[i])+1.7)
		g_array.append(0)
		b_array.append(0)

	fig=pylab.figure()
	pylab.xlabel("[Fe/H]")
        pylab.ylabel(r"$\rho/m_H$")
	#CP=pylab.loglog(dd,tk,'.',markersize=0.5,markerfacecolor=c_array,markeredgecolor=c_array)
	CP=pylab.scatter(r_array,dd,marker='o',s=3)
	pylab.yscale('log')
	pylab.axis([min(r_array),max(r_array),min(dd),max(dd)])
	pylab.savefig(filename)

	pylab.figure()
        pylab.ylabel("T/$\mu$")
        pylab.xlabel("[Fe/H]")
        pylab.scatter(r_array,tk,marker='o',s=3)
	pylab.yscale('log')
	pylab.axis([min(r_array),max(r_array),min(tk),max(tk)])
        pylab.savefig(filename2)
        
        fig=pylab.figure()
        pylab.ylabel("T/$\mu$")
        pylab.xlabel(r"$\rho/m_H$")
        CP=pylab.scatter(dd,tk,marker='o',s=3,c=r_array,alpha=0.5,edgecolors='none')
        CP_AX=pylab.clim([-5.,0.])
        pylab.colorbar(CP,CP_AX)
        pylab.xscale('log')
        pylab.yscale('log')
        pylab.axis([1e-1,1e6,1e1,1e4])
        pylab.savefig(filename3)


#!/usr/bin/env ipython3

##
# @file
# move centered positions to spherical coordinates

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

import numpy as np
import pdb
import gr_params as gpr

def run():
    for pop in range(gpr.pops):
        print('grh_Pos: input:')
        print(gpr.fileposcenter[pop])
        xall,yall,vlosall = np.loadtxt(gpr.fileposcenter[pop],
                                       comments='#', unpack=True)
        # 2*[Rscale], [km/s]

        Rall = np.sqrt(xall**2+yall**2) # [Rscale]
        rest = ( Rall <= gpr.Rmax ) # [Rscale]
        xall = xall[rest]           # [Rscale]
        yall = yall[rest]           # [Rscale]
        vlosall= vlosall[rest]      # [km/s]

        n = len(xall)
        x=xall                      # [Rscale]
        y=yall                      # [Rscale]
        vlos=vlosall                # [km/s]

        # old output to cartesian variables: needed anywhere?
        # print('grh_Pos: output:')
        # print(gpr.fileposcartesian)
        # fileposcartesian = open(gpr.fileposcartesian[pop], 'w')
        # print('# x [], y []',file=fileposcartesian)
        # for k in range(n):
        #     print(x[k], y[k], file=fileposcartesian)
        # fileposcartesian.close()

        # print(gpr.filevelcartesian)
        # filevelcartesian = open(gpr.filevelcartesian[pop], 'w')
        # print('# vlos [km/s]',file=filevelcartesian)
        # for k in range(n):
        #     print(vlos[k],file=filevelcartesian)
        # filevelcartesian.close()

        x = np.array(x);  y = np.array(y)   # 2*[Rscale]
        R = np.sqrt(x**2+y**2)              # [Rscale]

        frac = np.sum((R <= gpr.Rmax))*1./len(R)
        print('fraction of particles nearer than rmax: ',frac*100,'%')

        # phi (azimuthal angle [-pi,pi])
        Phi   = np.zeros(n)
        for k in range(n):
            if x[k]>0:                      # [Rscale]
                Phi[k] = np.arctan(y[k]/x[k])
            elif x[k]==0:                   # [Rscale]
                Phi[k] = np.sign(y[k])*np.pi/2
            elif y[k]>=0:                   # [Rscale]
                Phi[k] = np.arctan(y[k]/x[k])+np.pi
            elif y[k]<0:                    # [Rscale]
                Phi[k] = np.arctan(y[k]/x[k])-np.pi

        print(gpr.fileposspherical[pop])
        fileposspherical = open(gpr.fileposspherical[pop], 'w')
        print('# R [Rscale]', 'Phi [rad]', 'vlos [km/s]', file=fileposspherical)
        for k in range(n):
            print(R[k], Phi[k], vlos[k], file=fileposspherical)
        fileposspherical.close()
## \fn run()
# main functionality, to be called from gravlite

        
if __name__ == "__main__":
    gpr.showplots = True
    run()
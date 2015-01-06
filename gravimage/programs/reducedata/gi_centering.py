#!/usr/bin/env ipython3

##
# @file
# basic centering algorithms, shared between gr*_com

# (c) GPL v3 2014 ETHZ, Pascal S.P. Steger, psteger@phys.ethz.ch

import pdb
import numpy as np
import gi_helper as gh

def com_mean(x, y, z, pm):
    com_x = 1.*np.sum(x*pm)/np.sum(pm) # [pc]
    com_y = 1.*np.sum(y*pm)/np.sum(pm) # [pc]
    com_z = 1.*np.sum(z*pm)/np.sum(pm) # [pc]
    return com_x, com_y, com_z
## \fn com_mean(x,y,z,pm)
# mean COM, weighted by probability of membership
# @param x array of x values in [pc]
# @param y array of y values in [pc]
# @param z array of z values in [pc]
# @param pm array of probability of membership in [1]
# @return x,y,z of center of mass


def com_shrinkcircle_v_2D(x, y, vlos, pm):
    #eps = 1e-6
    com_x = 1.*np.sum(x*pm)/np.sum(pm);    com_y = 1.*np.sum(y*pm)/np.sum(pm);
    com_vlos = 1.*np.sum(vlos*pm)/np.sum(pm)
    bucom_x = com_x; bucom_y = com_y; bucom_vlos = com_vlos
    x -= com_x; y -= com_y; vlos -= com_vlos
    dr = np.sqrt(com_x**2+com_y**2)
    r0 = np.sqrt(x**2+y**2)

    nit = 0; minlen = len(x)/2.
    while nit < 200 and len(x) > minlen:
        nit += 1
        print('it ',nit,' with ',len(x),\
              ' part, COM=', \
              gh.pretty(bucom_x), gh.pretty(bucom_y),\
              ' offset ', gh.pretty(dr))

        # shrink sphere:
        # 1) calc radius
        r0 = np.sqrt(x**2+y**2)
        # 2) sort remaining particles
        order = np.argsort(r0)
        r0 = np.array(r0)[order]; x = np.array(x)[order]; y = np.array(y)[order]; pm = np.array(pm)[order]
        vlos = np.array(vlos)[order]

        # 3) cut x,y,z,pm after 1-10%
        end = len(r0)*0.95
        r0 = r0[:end]; x = x[:end]; y = y[:end]; vlos = vlos[:end]; pm = pm[:end]

        # calculate new COM
        com_x = 1.*np.sum(x*pm)/np.sum(pm);    com_y = 1.*np.sum(y*pm)/np.sum(pm)
        com_vlos = 1.*np.sum(vlos*pm)/np.sum(pm)
        dr = np.sqrt(com_x**2+com_y**2)

        # add to bucom
        bucom_x += com_x; bucom_y += com_y; bucom_vlos += com_vlos

        # recenter particles
        x -= com_x; y -= com_y; vlos -= com_vlos

    return bucom_x, bucom_y, bucom_vlos
## \fn com_shrinkcircle_v_2D(x, y, vlos)
# shrinking sphere in 2D, with LOS velocity
# @param x array of x values in [pc]
# @param y array of y values in [pc]
# @param vlos array of line of sight velocities in [km/s]



def com_shrinkcircle(x, y, z, pm):
    eps = 1e-6
    com_x = 1.*np.sum(x*pm)/np.sum(pm)
    com_y = 1.*np.sum(y*pm)/np.sum(pm)
    com_z = 1.*np.sum(z*pm)/np.sum(pm)
    bucom_x = 0.+com_x; bucom_y = 0.+com_y; bucom_z = 0.+com_z
    x -= com_x; y -= com_y; z -= com_z
    dr = np.sqrt(com_x**2+com_y**2+com_z**2)
    r0 = np.sqrt(x**2+y**2+z**2)

    nit = 0; minlen = len(x)*0.666666666
    while nit < 200 and len(x) > minlen:
        nit += 1
        print('it ',nit,' with ',len(x), ' part',\
              ' COM= ', gh.pretty(bucom_x), gh.pretty(bucom_y), gh.pretty(bucom_z),\
              ' offset ', gh.pretty(dr))

        # shrink sphere:
        # 1) calc radius
        r0 = np.sqrt(x**2+y**2+z**2)
        # 2) sort remaining particles
        order = np.argsort(r0)
        r0 = np.array(r0)[order]
        x = np.array(x)[order]
        y = np.array(y)[order]
        z = np.array(z)[order]
        pm = np.array(pm)[order]

        # 3) cut x,y,z,pm after 1-10%
        end = len(r0)*0.95
        r0 = r0[:end]; x = x[:end]; y = y[:end]; z = z[:end]; pm = pm[:end]

        # calculate new COM
        com_x = 1.*np.sum(x*pm)/np.sum(pm)
        com_y = 1.*np.sum(y*pm)/np.sum(pm)
        com_z = 1.*np.sum(z*pm)/np.sum(pm)
        dr = np.sqrt(com_x**2+com_y**2+com_z**2)

        # add to bucom
        bucom_x += com_x; bucom_y += com_y; bucom_z += com_z

        # recenter particles
        x -= com_x; y -= com_y; z -= com_z

    return bucom_x, bucom_y, bucom_z
## \fn com_shrinkcircle(x, y, z, pm)
# 3D shrinking sphere
# @param x array of x values in [pc]
# @param y array of y values in [pc]
# @param z array of z values in [pc]
# @param pm array of probability of membership values in [1]

def com_shrinkcircle_2D(x, y):
    com_x = np.mean(x)
    com_y = np.mean(y)
    bucom_x = 0.+com_x
    bucom_y = 0.+com_y
    x -= com_x
    y -= com_y
    dr = np.sqrt(com_x**2+com_y**2)
    R0 = np.sqrt(x**2+y**2)

    nit = 0; minlen = len(x)*0.666666666
    while nit < 200 and len(x) > minlen:
        nit += 1
        print('it ',nit,' with ',len(x), ' part',\
              ' COM= ', gh.pretty(bucom_x), gh.pretty(bucom_y),\
              ' offset ', gh.pretty(dr))

        # shrink sphere:
        # 1) calc radius
        R0 = np.sqrt(x**2+y**2)
        # 2) sort remaining particles
        order = np.argsort(R0)
        R0 = np.array(R0)[order]
        x = np.array(x)[order]
        y = np.array(y)[order]

        # 3) cut x,y,z,pm after 1-10%
        end = len(R0)*0.95
        R0 = R0[:end]
        x = x[:end]
        y = y[:end]

        # calculate new COM
        com_x = np.mean(x)
        com_y = np.mean(y)
        dr = np.sqrt(com_x**2+com_y**2)

        # add to bucom
        bucom_x += com_x
        bucom_y += com_y

        # recenter particles
        x -= com_x
        y -= com_y

    return bucom_x, bucom_y
## \fn com_shrinkcircle(x, y)
# 2D shrinking sphere
# @param x array of x values in [pc]
# @param y array of y values in [pc]

def com_shrinkcircle_v(x, y, z, vz, pm):
    eps = 1e-6
    com_x = 1.*np.sum(x*pm)/np.sum(pm)
    com_y = 1.*np.sum(y*pm)/np.sum(pm)
    com_z = 1.*np.sum(z*pm)/np.sum(pm)
    com_vz = 1.*np.sum(vz*pm)/np.sum(pm)
    bucom_x = 0.+com_x; bucom_y = 0.+com_y; bucom_z = 0.+com_z; bucom_vz = 0.+com_vz
    x -= com_x; y -= com_y; z -= com_z; vz -= com_vz
    dr = np.sqrt(com_x**2+com_y**2+com_z**2)
    r0 = np.sqrt(x**2+y**2+z**2)

    nit = 0; minlen = len(x)*0.666666666
    while nit < 200 and len(x) > minlen:
        nit += 1
        print('it ',nit,' with ',len(x), ' part',\
              ' COM= ', \
              gh.pretty(bucom_x), gh.pretty(bucom_y), gh.pretty(bucom_z),\
              ' vel=', gh.pretty(bucom_vz),\
              ' offset ', gh.pretty(dr))

        # shrink sphere:
        # 1) calc radius
        r0 = np.sqrt(x**2+y**2+z**2)
        # 2) sort remaining particles
        order = np.argsort(r0)
        r0 = np.array(r0)[order]
        x = np.array(x)[order]
        y = np.array(y)[order]
        z = np.array(z)[order]
        vz = np.array(vz)[order]
        pm = np.array(pm)[order]

        # 3) cut x,y,z,pm after 1-10%
        end = len(r0)*0.95
        r0 = r0[:end]; x = x[:end]; y = y[:end]; z = z[:end]; vz = vz[:end]; pm = pm[:end]

        # calculate new COM
        pmsum = np.sum(pm)
        com_x = 1.*np.sum(x*pm)/pmsum
        com_y = 1.*np.sum(y*pm)/pmsum
        com_z = 1.*np.sum(z*pm)/pmsum
        com_vz = 1.*np.sum(vz*pm)/pmsum
        dr = np.sqrt(com_x**2+com_y**2+com_z**2)

        # add to bucom
        bucom_x += com_x; bucom_y += com_y; bucom_z += com_z; bucom_vz += com_vz

        # recenter particles
        x -= com_x; y -= com_y; z -= com_z; vz -= com_vz

    return bucom_x, bucom_y, bucom_z, com_vz
## \fn com_shrinkcircle_v(x, y, z, vz, pm)
# shrinking sphere in 3D, following v_LOS
# @param x array of x values in [pc]
# @param y array of y values in [pc]
# @param z array of z values in [pc]
# @param vz array of LOS velocities (or rather, along z coordinate) in [km/s]
# @param pm array of probabilities of membership in [1]

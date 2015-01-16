#!/usr/bin/env ipython3

##
# @file
# unit definitions
# unit1__unit2 gives unit1 in terms of unit2
# naming convention: m^2/kg => m2kg_1

# (c) GPL v3 2014 Pascal Steger psteger@phys.ethz.ch

import numpy as np

day__hr = 24
hr__min = 60
day__min = day__hr*hr__min
min__s = 60
hr__s = hr__min*min__s
day__s = day__min*min__s
year__day = 365.25
year__s = year__day * day__s
century__s = 100 * year__s
arcsec__mas = 1000
full_circle__rad = 2*np.pi
full_circle__deg = 360
rad__deg = full_circle__deg/full_circle__rad
deg__arcmin = 60
rad__arcmin = rad__deg*deg__arcmin
arcmin__arcsec = 60
rad__arcsec = rad__arcmin*arcmin__arcsec
km__m = 1000


G1__m3kg_1s_2  = 6.67398e-11
pc__m  = 3.08567758e16
Msun__kg= 1.981e30
km__m  = 1000.
kpc__pc = 1000.
G1__pcMsun_1km2s_2  = G1__m3kg_1s_2*Msun__kg/km__m**2/pc__m

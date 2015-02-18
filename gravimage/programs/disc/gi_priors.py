#!/usr/bin/env ipython3

##
# @file
# check all parameters for prior constraints
# disc version

# (c) GPL v3 2015 ETHZ Pascal Steger, pascal@steger.aero

import pdb


def check_tilt(tmp_tilt, gp):
    return False
## \fn check_tilt(tmp_tilt, gp)
# check that tilt is physically possible
# @param tmp_tilt tilt profile
# @param gp


def check_bprior(rhocheck, nucheck):
    for jj in range(len(rhocheck)):
        if rhocheck[jj] < nucheck[jj]:
            return True
    return False
## \fn check_bprior(rhocheck, nucheck)
# check that observed tracer mass is less than total mass
# @param rhocheck density profile in [Munit/pc^3]
# @param nucheck density profile in [Munit/pc^3]

#!/usr/bin/env python3

##
# @file
# run all Walker data readout and binning routines in correct order

# (c) 2013 Pascal Steger, psteger@phys.ethz.ch

from grw_vardefs import *

import grw_com         # gives centerpos_?.txt
import gr_sphere
import grw_dmenclosed  # # gives dmenclosedmass_0.txt, totenclosedmass.txt
import grw_dens        # # gives densityfalloff_?.txt, enclosedmass_?.txt
import grw_siglos      # gives velocitydispersionlos_?.txt

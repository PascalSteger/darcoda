#!/usr/bin/env ipython-python3.2
import numpy as np
import numpy.random as npr
import pdb

import gl_params
import gl_plot as gpl
from gl_project import *
from gl_int import *
from gl_analytic import *

gpl.ion()
gpl.start()

r0 = np.arange(0.05,3.0,0.11)*1000
ascale = 1000.; Mscale =1.e6
nuanf = rho_anf(r0, ascale, Mscale)
rhoanf = rho_anf(r0, ascale, Mscale) * (1.+npr.uniform(-1.0,1.0,len(r0)))
gpl.plot(r0, rhoanf,'.',color='cyan', label=r'$\rho$')
gpl.plot(r0, rho_anf(r0,ascale,Mscale),color='cyan',lw=1)
betaanf = np.zeros(len(r0))
intbetaanf = np.zeros(len(r0))

# gpl.yscale('linear')

# surface density
Siganf = Sigma_anf(r0, ascale, Mscale)
gpl.plot(r0, Siganf, color='blue', lw=1)
Rhoself = rho_INT_Rho(r0, rhoanf)
gpl.plot(r0, Rhoself, '.', color='blue', label=r'$\Sigma$')
# gpl.plot(r0, M_anf(r0,1.,1.e6), color='blue')

# sigma_r^2
# gpl.plot(r0, sig_r_2_anf(r0, ascale, Mscale), color='green', lw=1)

siglosself, kappaself = ant_sigkaplos2surf(r0, betaanf, intbetaanf, rhoanf, nuanf)
gpl.plot(r0, siglosself, '.', color='orange', label=r'$\sigma_p*\Sigma$')
siganf = Sigma_sig_los_2_anf(r0, ascale, Mscale)
gpl.plot(r0, siganf, color='orange', lw=1)

# sigma_LOS^2
siglosanf = sig_los_anf(r0, ascale, Mscale)
gpl.plot(r0, siglosanf, color='black', lw=1)
gpl.plot(r0, np.sqrt(siglosself/Rhoself), '.', color='black', label=r'$\sigma_p$')

# Now add the legend with some customizations.
legend = gpl.legend(loc='lower left', shadow=True)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('0.90')

print('error fraction determined/true = ', (siglosself/Rhoself)/siglosanf)


# gpl.plot(r0, kappa_anf(r0))
# gpl.plot(r0, (kappaself/Rhoself)/(np.sqrt(siglosself/Rhoself)**4))

gpl.savefig('integration_siglos_disturbed.png')

gpl.ioff(); gpl.show()

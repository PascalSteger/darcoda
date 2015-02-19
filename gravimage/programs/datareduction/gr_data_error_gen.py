# Data_Error_Gen
# Adds Gaussian errors to (z, vz) data
# Method: take (z, vz) points for each star as centres of Gaussians with spread
#   sigma_z, sigma_vz. Then resample the star's (z, vz) data from that Gaussian

import numpy as np
import numpy.random as rand
import pdb

in_filename = '/home/hsilverw/LoDaM/darcoda/Data_Sets/simplenu/simplenu_sigz_raw.dat'
out_filename = '/home/hsilverw/LoDaM/darcoda/Data_Sets/simplenu/simplenu_sigz_raw_err_z_0p05_sdvz_5.dat'

error_z = 0.05 # fraction, eg 0.05 = 5%
sd_vx = 5.0 #[km s^-1]

data = np.loadtxt(in_filename)

z_data = data[:, 0] #[kpc]
vz_data = data[:, 1] #[km/s]

z_new_vec = []
vz_new_vec = []
new_data = []

for jter in range(0, len(z_data)):
    z_star = z_data[jter]
    vz_star = vz_data[jter]

    z_new = abs(rand.normal(loc=z_star, scale=z_star*error_z))
    vz_new = rand.normal(loc=vz_star, scale=sd_vx)

    z_new_vec.append(z_new)
    vz_new_vec.append(vz_new)
    new_data.append([z_new, vz_new])



z_new_vec = np.array(z_new_vec)
vz_new_vec = np.array(vz_new_vec)


np.savetxt(out_filename, new_data)

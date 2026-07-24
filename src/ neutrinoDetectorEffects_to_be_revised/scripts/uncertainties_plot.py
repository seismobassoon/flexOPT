import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gmean

import sys
import os
path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(path, 'scripts'))
from plot_tools import *

path_project = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
folder_data  = os.path.join(path_project, 'data/')
file    = 'DUNE-prem64-OutputExpTruth.root'
channel = 'h_tracks_all'
histo   = SyntheticData(folder_data, file, channel)
x, y, z = histo.get_histo()
ct_bins, enu_bins = z.shape

def sigma_E_rel(Enu, Eres_cst, Eres_sqrt):
    return Eres_cst + Eres_sqrt / np.sqrt(Enu)

def sigma_E(Enu, Eres_cst, Eres_sqrt):
    return (Eres_cst + Eres_sqrt / np.sqrt(Enu)) * Enu

def sigma_ct(ct, ctres_cst, ctres_sqrt):
    return ctres_cst + ctres_sqrt / np.sqrt(np.arccos(ct))

# ORCA
Eres_cst   = 0.25
Eres_sqrt  = 0.0
ctres_cst  = 0.0 * (np.pi/180)
ctres_sqrt = 30.0 * (np.pi/180)
# HK
Eres_cst   = 0.15
Eres_sqrt  = 0.0
ctres_cst  = 0.0 * (np.pi/180)
ctres_sqrt = 15.0 * (np.pi/180)
# DUNE
Eres_cst   = 0.05
Eres_sqrt  = 0.0
ctres_cst  = 5.0 * (np.pi/180)
ctres_sqrt = 0.0 * (np.pi/180)


var_z = np.zeros((ct_bins, enu_bins))
for i in range(ct_bins):
    for j in range(enu_bins):
        E_mean = gmean([x[j], x[j+1]])
        ct_mean = np.mean([y[i], y[i+1]])
        Eres  = sigma_E(E_mean, Eres_cst, Eres_sqrt)
        ctres = sigma_ct(ct_mean, ctres_cst, ctres_sqrt)
        var_z[i, j] = Eres**2 + ctres**2

plt.rcParams.update({'font.size': 18})
fig, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
im      = ax.pcolormesh(x, y, var_z)
ax.set_xlabel(r"$E_\nu$ [GeV]")
ax.set_ylabel(r"$\cos \theta_z$")
ax.set_xscale('log')
ax.set_xticks([1, 2, 4, 6, 8, 10, 20, 30, 40])
ax.set_xticklabels([1, 2, 4, 6, 8, 10, 20, 30, 40])
cbar = fig.colorbar(im, ax=ax, orientation='horizontal', location='top')
cbar.set_label(r"$\sigma_E^2 + \sigma_\theta^2$ (E in GeV, $\theta$ in radians)", labelpad=10)
plt.show()



#RMS(E)/E    = E_res_cst + E_res_sqrt / sqrt(E/GeV)
#RMS(zenith) = Theta_res_cst + Theta_res_sqrt / sqrt(E/GeV)



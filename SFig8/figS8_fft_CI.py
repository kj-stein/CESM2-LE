#README
#This is shared FOR REFERENCE ONLY, as supporting material for the paper 
#"Ubiquity of human-induced changes in climate variability", 
#Rodger et al, Earth System Dynamics, 2021.
#(https://doi.org/10.5194/esd-2021-50)
#
#This shows how the confidence interval in FigS8 is generated. 
#
#Authors: Ji-Eun Kim (jieunkim[at]pusan[dot]ac[dot]kr)


import numpy as np
import xarray as xr 
import matplotlib.pyplot as plt 


def plot_fft_ci():
    #Plot FFT spectrum from a netcdf file (spectrum as function of [freq,lat,lon,ensemble]).

    dirsave = '/MY_DIRECTORY/'
    stryy = str(year1)+'-'+str(year2)
    spec_qyxe = xr.open_dataset(dirsave+'PRECT.spec_qyxe.nc')

    frq = spec_qyxe.frq.values
    prd = 1/frq

    weights = np.cos(np.radians(spec_qyxe.lat)) ; weights.name = 'weights'
    specW_qyxe = spec_qyxe.weighted(weights)
    spec_q = specW_qyxe.mean(['lat','lon','ensA'])[var].values

    spec_qe = specW_qyxe.mean(['lat','lon'])[var]

    #To get confidence interval, standard_error * z
    #sample size = 100 (100 ensembles)
    ci_q = spec_qe.std('ensA').values / 10 * 1.96 #95% confidence level

    #nlon = len(spec_qyxe.lon.values)
    #nlat = len(spec_qyxe.lat.values)
    #ci_q2 = spec_qyxe.std(['lat','lon','ensA'])[var].values /10/np.sqrt(nlat*nlon)*1.96


    #==========================================================FIG
    fig, ax = plt.subplots(figsize=(5,4))
    ax.plot(prd,spec_q, color='r', linewidth=1)
    ax.fill_between(prd, spec_q+ci_q, spec_q-ci_q, alpha=0.3,color='gray')

    plt.xscale('log')
    plt.yscale('log')

    plt.xlim(max(prd[1:]),2)
    plt.ylim(3e-2, 6e0)

    plt.xlabel('Period (day)', fontsize=20)
    plt.ylabel('Amp. (mm/day)', fontsize=20)

    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)

    plt.tight_layout()
    plt.show()

    fig.savefig('fft_ci.png', transparent=True, dpi=300)
    plt.close(fig)

    return None

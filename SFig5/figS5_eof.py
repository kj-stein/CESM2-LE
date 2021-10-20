#README
#This is shared FOR REFERENCE ONLY, as supporting material for the paper 
#"Ubiquity of human-induced changes in climate variability", 
#Rodger et al, Earth System Dynamics, 2021.
#(https://doi.org/10.5194/esd-2021-50)
#
#This calculates PCs and EOFs and plot them.
#
#Authors: Ji-Eun Kim (jieunkim[at]pusan[dot]ac[dot]kr)


import numpy as np
import xarray as xr 
import matplotlib.pyplot as plt 
from scipy import signal

#Please see Fig1 codes for FigS5a.

class EOFplot_yx:
    def __init__(self, data_tyx,time,lat,lon,npc=1,detrend=False):
        if detrend: data_tyx = signal.detrend(data_tyx,axis=0)
        self.data_tyx = data_tyx
        self.time = time
        self.lat = lat
        self.lon = lon
        self.npc = npc
 

    def eofsolver(self):
        weights_array = np.sqrt(np.cos(np.deg2rad(self.lat)))[:, np.newaxis]
        solver = MultivariateEof([self.data_tyx], weights=[weights_array])
        vfrac = solver.varianceFraction(neigs=self.npc)*100
    
        pcs = solver.pcs(npcs=self.npc, pcscaling=1)
        eofs = np.asarray(solver.eofs(neofs=self.npc, eofscaling=2))

        self.vfrac = vfrac
        self.pcs = pcs
        self.eofs = eofs

        return #vfrac,pcs,eofs


    def eofplot(self):
        for pii in range(self.npc):
            eof_yx = self.eofs[0,pii,:,:]
            pc_tt = self.pcs[:,pii] 
            
            spii = str(pii+1)
            figname = self.'PC'+spii+'.png'
            figtit1 = self.'PC'+spii+' {:2.0f}'.format(self.vfrac[pii])+'%'
            figtit2 = self.'EOF'+spii+' {:2.0f}'.format(self.vfrac[pii])+'%'

            projection=ccrs.PlateCarree(central_longitude=180)
            transform=ccrs.PlateCarree()
        
            fig, (ax1,ax2) = plt.subplots(2,1,figsize=(7,7))
        
            ax1 = plt.subplot(211)
            ax1.plot(time,pc_tt)
            ax1.set_ylabel('PC')
            ax1.set_title(figtit1)
        
            ax2 = plt.subplot(212,projection=projection)
            cf = ax2.contourf(lon,lat,eof_yx,transform=transform)
            plt.colorbar(cf, orientation='vertical',fraction=0.013,aspect=35,pad=0.08)
            ax2.coastlines()
            ax2.set_title(figtit2)
            
            fig.savefig(figname, transparent=False,dpi=300)
            #plt.show()
            plt.close(fig)

        return None

#Example with data array of data3D_tyx,time,lat,lon
#obj = EOFplot_yx(data3D_tyx,time,lat,lon,npc=1,detrend=True)
#obj.eofsolver()
#obj.eofplot()

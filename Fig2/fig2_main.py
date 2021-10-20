#README
#This is shared FOR REFERENCE ONLY, as supporting material for the paper 
#"Ubiquity of human-induced changes in climate variability", 
#Rodger et al, Earth System Dynamics, 2021.
#(https://doi.org/10.5194/esd-2021-50)
#
#This python script shows how Figure 2 and Figure S7 (map with spectrum and PDF) are generated.
#The code requires pre-processed data netcdf and csv files.
#
#Authors: Ji-Eun Kim (jieunkim[at]pusan[dot]ac[dot]kr) for map and spectrum
#         Lei Hwang  (huanglei[at]pusan[dot]ac[dot]kr) for PDF


import numpy as np
import pandas as pd
import xarray as xr 
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
import cartopy.crs as ccrs
import cmocean


def simpledrawbox(pt1,pt2,color='black',linewidth=2,shiftx=0):
   #Draws an area box.
   transform=ccrs.PlateCarree()
  
   xx1 = pt1[0] ; yy1 = pt1[1]
   xx2 = pt2[0] ; yy2 = pt2[1]
 
   xx1 = xx1-shiftx
   xx2 = xx2-shiftx
 
   plt.plot([xx1,xx2],[yy1,yy1], color=color, linewidth=linewidth,transform=transform)
   plt.plot([xx1,xx2],[yy2,yy2], color=color, linewidth=linewidth,transform=transform)
   plt.plot([xx1,xx1],[yy1,yy2], color=color, linewidth=linewidth,transform=transform)
   plt.plot([xx2,xx2],[yy1,yy2], color=color, linewidth=linewidth,transform=transform)
   return None


def plot_map(ax):
    #Plots the central map of surface temperature and precipitation difference.
    dr = '/MY_DIRECTORY/'
    #Open netcdf files of difference data (2070-2099)-(1960-1989) as xarray Dataset.
    yx1 = xr.open_dataset(dr+'TS.diff.ens100.ANN.nc')   #surface temperature 
    yx2 = xr.open_dataset(dr+'PRECT.diff.ens100.ANN.nc')#total precipitation
    lon = yx1.lon
    lat = yx2.lat
    yx1 = yx1.TS
    yx2 = yx2.PRECT

    lon = np.concatenate((lon,np.array([360]))) #for continuity at 0 degree longitude.
    yx1 = np.concatenate((yx1,yx1[:,:1]),axis=1)
    yx2 = np.concatenate((yx2,yx2[:,:1]),axis=1)

    transform=ccrs.PlateCarree()

    lev1 = [0,1,2,3,4,5,6,7,8,9,10,11] ; cmap1 = plt.cm.hot.reversed()
    cf = ax.contourf(lon,lat,yx1,levels=lev1,cmap=cmap1,transform=transform,extend='both')

    grid_x, grid_y = np.meshgrid(lon[0::5], lat[0::5]) #for dots
    x, y = grid_x.flatten(), grid_y.flatten()
    yx2 = yx2[0::5,0::5]
    #for positive precipitation diff
    tmp = np.copy(yx2) 
    tmp[np.where(tmp < 0)] = 0
    area = tmp.flatten()*8.5
    sct1 = ax.scatter(x,y,s=area, c='royalblue',transform=transform)
    handles, labels = sct1.legend_elements(num=[6,12,18,24], prop="sizes", color='royalblue')#, alpha=0.6)
    ax.legend(handles, labels, loc='lower right', bbox_to_anchor=(1.04,0.0), title="")
    legend1 = ax.legend(handles, labels, loc= (0.93,0.73),  
                    title='wet',
                    labelspacing = 0.1, 
                    handletextpad = 0.2,
                    frameon = False
                    )   
    ax.add_artist(legend1)

    #for negative precipitation diff
    yx2[np.where(yx2 > 0)] = np.nan
    area = -yx2.flatten()*8.5
    sct2 = ax.scatter(x,y,s=area, c='cyan',transform=transform)
    kw = dict(prop = 'sizes', num = (6,12,18,24), color = 'cyan')
    ax.legend(*sct1.legend_elements(**kw), #use the same marker size for 'wet' and 'dry'
                    loc = (0.93,0.0),
                    title = 'dry',
                    labelspacing = 0.1,
                    handletextpad = 0.2,
                    frameon = False
                    )

    ax.coastlines()
    plt.colorbar(cf, orientation='vertical',fraction=0.015,aspect=25,pad=0.06)

    simpledrawbox([190,-5],[240,5],color='k')   #nino
    simpledrawbox([235,32],[242,41],color='k')  #california
    simpledrawbox([300,40],[345,60],color='k')  #north atlantic
    simpledrawbox([280,-10],[310,10],color='k') #amazon

    return None


def rd_spec(var,freq,ensGR,region,year1,year2,anm='raw'):
    #Reads a netcdf spectrum file of a spectrum.
    dr = '/MY_DIRECTORY/'
    stryy = str(year1)+'-'+str(year2)
    spec_q = xr.open_dataset(dr+'Spectrum_amp.'+var+'.'+region+'.'+ensGR+'.'+stryy+'.nc')
    return spec_q[var].values,spec_q.frq.values




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#From here, it's the main.-----------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Read spectrum data as xarray Dataset.
nfire1,frq = rd_spec('NFIRE','day_1','ens90','California_land',1960,1989)
nfire2,frq = rd_spec('NFIRE','day_1','ens90','California_land',2070,2099)
chl1,frq = rd_spec('allChl_SURF','day_1','ens100','NAtlantic',1960,1989)
chl2,frq = rd_spec('allChl_SURF','day_1','ens100','NAtlantic',2070,2099)
prect1,frq = rd_spec('PRECT','day_1','ens100','Nino3.4',1960,1989)
prect2,frq = rd_spec('PRECT','day_1','ens100','Nino3.4',2070,2099)
sst1,frq = rd_spec('SST','day_1','ens100','Nino3.4',1960,1989)
sst2,frq = rd_spec('SST','day_1','ens100','Nino3.4',2070,2099)
nep1,frq = rd_spec('NEP','day_1','ens90','Amazon_land',1960,1989)
nep2,frq = rd_spec('NEP','day_1','ens90','Amazon_land',2070,2099)

#Read PDF (histogram) data (saved as csv).
calif_nfire_raw  = pd.read_csv('/proj/lhuang/LENS/Keith_fig2_PDFs/10year_100ens_pdf/nfire_calif_raw_PDFs.csv')
na_chl_raw       = pd.read_csv('/proj/lhuang/LENS/Keith_fig2_PDFs/10year_100ens_pdf/chl_NA_raw_PDFs.csv')
prect_nino_raw   =pd.read_csv('/proj/lhuang/LENS/Keith_fig2_PDFs/10year_100ens_pdf/prect_nino_raw_PDFs.csv')
sst_nino_raw     = pd.read_csv('/proj/lhuang/LENS/Keith_fig2_PDFs/10year_100ens_pdf/SST_nino3_raw_PDFs.csv')
nep_amazon_raw   = pd.read_csv('/proj/lhuang/LENS/Keith_fig2_PDFs/10year_100ens_pdf/nep_amazon_raw_PDFs.csv')

data_sp = [calif_nfire_raw,na_chl_raw,nep_amazon_raw,prect_nino_raw,sst_nino_raw]
xlabel_sp = ['FireCounts California_land\n[counts/(100km)$^2$/yr]',
           'Chlorophyll NAtlantic\n(mg/m$^3$)',
           'NEP Amazon_land\n(gC/km$^2$/s)',
           'Precipitation Nino3.4\n(mm/day)',
           'SST Nino3.4\n(\u2103)']
xlabel_sp = ['Fire (counts/(100km)$^2$/yr)',
           'Chlorophyll (mg/m$^3$)',
           'NEP (gC/km$^2$/s)',
           'Precipitation (mm/day)',
           'SST (\u2103)']

plt.rcParams.update({'font.size':12.0})
color_hist, color_proj = 'C9', 'C3'

ylogscale = ['log','log','linear','log','linear']
minor_bool = [True,True,False,True,False]
alphas = [r'$\mu:\ $',r'$\sigma^{2}:\ $', r'$\tilde{\mu}_{3}:\ $', r'$\kappa:\ $']
stats_inx = [4,8,16,6,10]
anno_hist_x = [0.03,0.28,0.08,0.15,0.05]
anno_proj_x = [0.54,0.63,0.58,0.58,0.05]
anno_y = [0.9,0.9,0.9,0.9,0.5]
xlim_l = [-50,-0.5,-50,-15,16]
xlim_r = [2200,18,50,800,36]

fig, ax = plt.subplots(1,1,
                       figsize = (18,14),
                       subplot_kw={'projection':ccrs.Robinson(central_longitude=180)})

plt.subplots_adjust(left = 0.29, bottom = 0.297, right = 0.69, top = 0.697)
plot_map(ax)

#Add axes for subplots of spectrum and PDF.
ax1 = fig.add_axes([0.14,0.69,0.14,0.14])
ax2 = fig.add_axes([0.14+0.14+0.05,0.69,0.14,0.14])
ax3 = fig.add_axes([0.15+0.14+0.14+0.05+0.05+0.02,0.69,0.14,0.14])
ax4 = fig.add_axes([0.15+0.14+0.14+0.14+0.05+0.05+0.05+0.02,0.69,0.14,0.14])
ax5 = fig.add_axes([0.76, 0.28 + 0.14 + 0.02, 0.14,0.14])
ax6 = fig.add_axes([0.76, 0.28 - 0.04, 0.14,0.14])
ax7 = fig.add_axes([0.73 - 0.14-0.05 - 0.14-0.05,0.18, 0.14,0.14])
ax8 = fig.add_axes([0.73 - 0.14-0.05,0.18, 0.14,0.14])
ax9  = fig.add_axes([0.12, 0.28 + 0.14 + 0.02,0.14,0.14])
ax10 = fig.add_axes([0.12, 0.28 - 0.04,0.14,0.14])


pdf_ax_list = [ax2,ax4,ax6,ax8,ax10]
for i in range(5):
    ax = pdf_ax_list[i]
    ax.set_ylabel('Probability Density', labelpad = 0.01)
    ax.set_xlabel(xlabel_sp[i], labelpad = 0.01, linespacing = 1)

    ax.plot(data_sp[i]['bins'], data_sp[i]['h_hist'], 'o', markersize=1.5, color = color_hist)
    ax.plot(data_sp[i]['bins'], data_sp[i]['h_proj'], 'o', markersize=1.5, color = color_proj)

    ax.set_xlim(xlim_l[i],xlim_r[i])
    ax.set_yscale(ylogscale[i])
    ax.minorticks_on()
    if minor_bool[i]:
        y_minor = mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1, 10) * 0.1, numticks = 10)
        ax.yaxis.set_minor_locator(y_minor)
        ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

spec_ax_list = [ax1,ax3,ax5,ax7,ax9]
prd = 1/frq
specs1 = [nfire1, chl1, nep1 ,prect1, sst1]
specs2 = [nfire2, chl2, nep2 ,prect2, sst2]
ylim_l = [9e-2  , 2e-4, 1e-1, 3e-2  , 1e-4]
ylim_r = [2e2   , 2e0 , 1e1 , 6e0   , 3e0 ]

titles = ['California Fire Counts',
          'N. Atlantic Chlorophyll',
          'Amazon NEP',
          'Nino3.4 Precipitation',
          'Nino3.4 SST']

for i in range(5):
    ax = spec_ax_list[i]
    ax.set_title(titles[i])

    ax.set_xlabel('Period (day)', linespacing = 1)
    ax.set_ylabel('Spectral Amplitude', labelpad = 0.01)

    ax.plot(prd,specs2[i], color=color_proj,linewidth=1.8)
    ax.plot(prd,specs1[i], color=color_hist,linewidth=1.5)
    if i == 1: ax.plot(prd,specs2[i], color=color_proj,linewidth=1.8)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlim(365*30, 2)
    ax.set_ylim(ylim_l[i],ylim_r[i])


#Draw grey boxes and lines. 
ax0 = fig.add_axes([0,0,1,1])
ax0.set_frame_on(False)
ax0.axes.get_xaxis().set_visible(False)
ax0.axes.get_yaxis().set_visible(False)

x1,x2,y1,y2 = 9.5,48.5,64.5,86    ; poly1 = [[x1,y1],[x2,y1],[x2,y2],[x1,y2]]
x1,x2,y1,y2 = 50.1,89.4,64.5,86   ; poly2 = [[x1,y1],[x2,y1],[x2,y2],[x1,y2]]
x1,x2,y1,y2 = 71.4,91.2,19.6,61.3 ; poly3 = [[x1,y1],[x2,y1],[x2,y2],[x1,y2]]
x1,x2,y1,y2 = 30.5,69.5,13.4,34.9 ; poly4 = [[x1,y1],[x2,y1],[x2,y2],[x1,y2]]
x1,x2,y1,y2 = 7.5,27.3,19.7,61.3  ; poly5 = [[x1,y1],[x2,y1],[x2,y2],[x1,y2]]

ax0.add_patch(Polygon(poly1, facecolor='None', edgecolor='grey', linewidth=2.5, alpha=0.6))
ax0.add_patch(Polygon(poly2, facecolor='None', edgecolor='grey', linewidth=2.5, alpha=0.6))
ax0.add_patch(Polygon(poly3, facecolor='None', edgecolor='grey', linewidth=2.5, alpha=0.6))
ax0.add_patch(Polygon(poly4, facecolor='None', edgecolor='grey', linewidth=2.5, alpha=0.6))
ax0.add_patch(Polygon(poly5, facecolor='None', edgecolor='grey', linewidth=2.5, alpha=0.6))

x1,x2,y1,y2 = 52.6,46,56,64.4
ax0.add_line(Line2D([x1,x2],[y1,y2], color='grey', linewidth=2.5, alpha=0.6))
x1,x2,y1,y2 = 60.3,62.5,58.8,64.4
ax0.add_line(Line2D([x1,x2],[y1,y2], color='grey', linewidth=2.5, alpha=0.6))
x1,x2,y1,y2 = 61,71.3,49.3,49.3
ax0.add_line(Line2D([x1,x2],[y1,y2], color='grey', linewidth=2.5, alpha=0.6))
x1,x2,y1,y2 = 51,51,35,48.9
ax0.add_line(Line2D([x1,x2],[y1,y2], color='grey', linewidth=2.5, alpha=0.6))
x1,x2,y1,y2 = 27.4,30.4,30,30
ax0.add_line(Line2D([x1,x2],[y1,y2], color='grey', linewidth=2.5, alpha=0.6))

ax0.text(63.3,37.0,'(mm/day) (\u2103)', fontsize=11.5)

ax0.set_xlim([0,100])
ax0.set_ylim([0,100])

#plt.show()
fig.savefig('fig2.pdf', transparent=False,dpi=400)

#README
#This is shared FOR REFERENCE ONLY, as supporting material for the paper 
#"Ubiquity of human-induced changes in climate variability", 
#Rodger et al, Earth System Dynamics, 2021.
#(https://doi.org/10.5194/esd-2021-50)
#
#This python script caculates an amplitude spectrum.
#
#Authors: Ji-Eun Kim (jieunkim[at]pusan[dot]ac[dot]kr)


import numpy as np
import matplotlib.pyplot as plt 


def fft1d(time,sig,ffttype='amp'):
    #Calculate Amplitude Spectrum from 1D array.
    #Output unit is the same as input signal unit.

    n = len(time) # length of the signal
    dt=abs(time[0]-time[1])

    Spec = np.fft.fft(sig)/n # fft computing and normalization
    frq = np.fft.fftfreq(n, d=dt) #two-side freq

    spec = 2*abs(Spec[:int(n/2)]) #take the positive frequency side only and double it.
    spec[0] =spec[0]/2 #No doubling for zero frequency
    frq = frq[:int(n/2)] #one-side freq

    return frq,spec

#e.g.
#time = np.arange(100) ; sig = 100*np.sin(2*np.pi/20 * time)
#frq,spec = fft1d(time,sig)
#fig, (ax1,ax2) = plt.subplots(2)
#ax1.plot(time,sig)
#ax2.plot(frq,spec)
#plt.show()

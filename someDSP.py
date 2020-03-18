#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 11:13:47 2020

@author: bazzi
"""

##############################################################################
##############################################################################

def getCWT(y,Fs,Fmin = 2 , Fmax = 140, res = 140, dec = 10, filt = 1, view = 0):
    
    import numpy as np
    import scipy as sc
    from scipy import signal
    
    # zero mean
    y = y - np.mean(y)
           
    # -optional- decimate (unless input dec = 0, in which case there is no decimation)
    
    if dec > 1:       
        y = sc.signal.decimate(y, dec, n=None, ftype='iir', axis=-1, zero_phase=True)
        Fs = Fs/dec
        
    # -optional- IIR fliter (butterworth, 120Hz lowpass)
        
    if filt == 1:
        
        # filter design
        b, a = signal.butter(10, 120, 'low',  analog=False, output='ba', fs=Fs)
        
        # apply filter
        #zi = signal.lfilter_zi(b, a)
        #z, _ = signal.lfilter(b, a, tlf2, zi=zi*tlf2[0])
        #z2, _ = signal.lfilter(b, a, z, zi=zi*z[0])

        y = signal.filtfilt(b, a, y)
    
    # set the frequency array with linear intervals (or log, commented out)

    # freq = np.space(np.log10(Fmin),np.log10(Fmax), num = res)     
    if res > Fmax:
        res = Fmax       
    freq = np.linspace(Fmin,Fmax,num=res-1,endpoint=True,dtype=int)
    
    # compute the CWTM   
    w = 8
    widths = w*Fs / (2*freq*np.pi)
    cwtm = signal.cwt(y, signal.morlet2, widths, w=w)
    
    # optional plot
    if view == 1:
        import seaborn as sns
        import pandas as pd 
        import matplotlib.pyplot as plt 
            
            
        tl = np.linspace(0, len(y)/Fs, len(y))
        df = pd.DataFrame(data=cwtm)
        df.columns = tl
        df.rows = freq
         
        figure, axes = plt.subplots(nrows=2, ncols=1) # sharex=True
        
        axes[0].set_title('raw data/cMor coeffs')
        axes[0].set(xlabel='', ylabel='mV')
        axes[0].plot(tl,y) # other way: ax1 = axes[0].plot(tl,tlf)
        
        #axes[1].set_title('cMor coeffs')
        axes[1] = sns.heatmap(abs(df),cbar=False,xticklabels=False,linewidths=0)
        axes[1].set(xlabel='time(s)', ylabel='freqs(Hz)')
        plt.show()                   
        return cwtm, freq, y, Fs, figure    

    else:
        return cwtm, freq, y, Fs
        
##############################################################################
##############################################################################
        
def getMI(y,Fs, res = 140, dec = 10, filt = 1, nPerm = 0, nJobs = 0):
    import numpy as np
    
    # get the CWT of the signal from the function above
    cwtm, freq, y, Fs = getCWT(y,Fs,2,140,res,dec,1,0)
    
    # find the indexes of the closest values to 2 and 30 (Hz) in the freq array 
    # (necessary for log freq array and different freq resolutions)
    Lfreq = np.zeros(2)
    val, Lfreq[0] = min((val, idx) for (idx, val) in enumerate(freq - 2)) 
    val, Lfreq[1] = min((val, idx) for (idx, val) in enumerate(abs(30 - freq)))
    
    # find the indexes of the closest values to 30 and 120 (Hz) in the freq array
    Hfreq = np.zeros(2)
    val, Hfreq[0] = min((val, idx) for (idx, val) in enumerate(abs(30 - freq)))
    val, Hfreq[1] = min((val, idx) for (idx, val) in enumerate(120 - freq))
    
    # array of indexes of the Low and High freq based on the indexes above (necessary for the log option)
    phM = np.linspace(int(Lfreq[0]),int(Lfreq[1]),int(Lfreq[1]-Lfreq[0]), endpoint=True, dtype=int)
    amM = np.linspace(int(Hfreq[0]),int(Hfreq[1]),int(Hfreq[1]-Hfreq[0]), endpoint=True, dtype=int)
    
    # get the actual frequencies for low/high bands for export
    LfreqE = freq[int(Lfreq[0]) : int(Lfreq[1])]
    HfreqE = freq[int(Hfreq[0]) : int(Hfreq[1])]
    
    # calculate the phase and amplitude-envelope of the low and high freq bands 
    phi = np.arctan2(np.imag(cwtm[phM]),np.real(cwtm[phM])) # phase of the LF
    amp = abs(np.real(cwtm[amM])) # amplitude enevelope HF
        
    # preallocate
    CWTampSH = np.nan * np.ones(shape=(nPerm,len(HfreqE),np.size(cwtm,1)))
    MIsh = np.nan * np.ones(shape=(nPerm,len(HfreqE),len(LfreqE)))
    #tl = np.linspace(0,np.size(amp,1)/Fs,np.size(amp,1))
   
    # create surrogate dataset with shuffled amplitudes and calculate MI         
    for mm in range(nPerm): # phase array        
        for kk in range(len(amM)):
            CWTampSH[mm,kk,:] = np.random.permutation(amp[mm]) #shuffle the amplitude
            
         
    # set-up and execute MIcalc in multi processing (with the timer)
    # https://buildmedia.readthedocs.org/media/pdf/joblib/latest/joblib.pdf
    from joblib import Parallel, delayed, parallel_backend
    from joblib.externals.loky import set_loky_pickler
    import time
 
    start = time.time()
    
    with parallel_backend("multiprocessing"):
        set_loky_pickler('cloudpickle')
        MIsh = Parallel(n_jobs = nJobs)(
        delayed(calcMI)(phi,CWTampSH[nn],phM,amM)
        for nn in range(nPerm))
    
    stop = time.time()
    print('MI calculation on shuffled dataset: {:.2f} s'
      .format(stop - start))
          
    # create threshold matrix and kill non-significant coupling coefficients
    MImean = np.mean(MIsh,0);
    MIsd = np.std(MIsh,0); 
    
    import scipy.stats as sct
    
    thrMX = sct.norm.ppf(1-0.001,MImean,MIsd) # Threshold matrix: p = 0.001 inverse t distr
    thrMX = np.nan_to_num(thrMX)
    
    ModIdx = calcMI(phi,amp,phM,amM) # MI of the experimental dataset

    MI = ((ModIdx*(ModIdx > thrMX))) #MI of the dataset thresholded at 1% with the bootstrapped dataset
    
    
    return MI, LfreqE, HfreqE, y, Fs

##############################################################################
##############################################################################

    
def calcMI(phi,amp,phM,amM):   
    
    import numpy as np
        
    # preallocate
    ModIndx = np.nan * np.ones(shape=(np.size(amp,0),np.size(phi,0))) # preallocate MI matrix
    
    # numebr of bins for the phase
    N = 18 
       
    for jj in np.arange(0,len(phM)):
    
        for kk in np.arange(1,len(amM)): 
            
            # TortMethod
            
            phase = phi[jj,:] + np.pi
            amplitude = amp[kk,:]
                    
            rng = 360/N;     
            Afa = np.zeros((N,1))      
            
            for xx in np.arange(1,N):
                ST = (np.pi/180)*((xx-1)*rng)
                ET = (np.pi/180)*(xx)*rng      
                states = np.logical_and(phase>ST,phase<=ET) # stats contains a boolean that is true when both conditions aare true
                ind = np.where(states)[0] # get the indexes of the True(s) in states                        
                Afa[xx] = np.nanmean(amplitude[ind])
                
                
            # Calculate entropy measure H
            Afa = Afa[np.nonzero(Afa)]
            p = Afa/np.nansum(Afa)
            H = -np.nansum(p*np.log(p))
    
            # Find MI
            MI = (np.log(N) - H)/np.log(N);
            ModIndx[kk,jj]=MI
            
        
        
           
    return ModIndx   

        
##############################################################################    
##############################################################################
 
    
def showMI(MI, LfreqE, HfreqE, y, Fs):
       
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    
    tl = np.linspace(0, len(y)/Fs, len(y))
    df = pd.DataFrame(data = MI, index = HfreqE, columns = LfreqE)
        
    figure, axes = plt.subplots(nrows=2, ncols=1) 
        
    axes[0].set_title('raw data')
    axes[0].set(xlabel='', ylabel='mV')
    axes[0].plot(tl,y) # other way: ax1 = axes[0].plot(tl,tlf)
    
    #axes[1].set_title('cMor coeffs')
    axes[1] = sns.heatmap(abs(df),cbar=False)
    axes[1].set(xlabel='freqs(Hz)', ylabel='freqs(Hz)')
               
    
    
   
##############################################################################    
##############################################################################
        

    
    
    
    
    
    
    
    
    
        
        
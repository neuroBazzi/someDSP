# someDSP
some digital signla processing routines for neuroscience (CWT, Phase amplitude coupling)
someDSP intructions

Paolo Bazzigaluppi, paolo.bazzigaluppi@gmail.com

Some routines for digital signal processing (DSP) for neuroscience (e.g. LFP, EEG). At this time, it contains routines to estimate the Continous Wavelet Transform (cMorl wavelet) and Phase Amplitude Couling (i.e.  Modulation Index as in Tort et al. 2010.)

Use: 

from someDSP import getCWT, getMI, showMI

#################################################################################
- getGWT:  cwtm, freq, y, Fs = getCWT(y,Fs,Fmin,Fmax,res,dec,filt,view)

Inputs: y (your signal, as array), Fs (sampling frequency, integer), Fmin (lower bound for the CWT computation, integer, default = 2), Fmax (higher bound for the CWT computation, integer, default = 140), res (resolution of the CWT, integer, default = 140, cofficients for the cwt are calulated linearly between Fmin and Fmax, logistic distribution is avialable too, just comment/uncomment lines 43/44), dec (decimate factor if you want to downsample, integer, default = 10, set as 0 if you do not want to decimate), filt (digital filter, boolean, 1 = yes, 0 = no, default = 1, this will apply a 10 poles lowpass Butterworth IIR forward-backward filter with cutoff at 120Hz, you can change parameters at line 32), view (plotting, boolean, 1 for plot, 0 for no plot, default = 0). The outputs are the CWT coefficients, the frequencies array and the raw data with sampling frequency (these are useful in case of downlasmpling/filtering).

Comments: this function will be automatically called by getMI, if you use getMI you do not need to 
call this. 

#################################################################################
- getMI:  MI, LfreqE, HfreqE, y, Fs = (y,Fs, res = 140, dec = 10, filt = 1, nPerm = 0, nJobs = 0)

This rountine uses the getCWT to compute the complex wavelet transform of the input signal, and uses arctan of the real part of the transform and its complex conjugate to estimate the istantenous phase. It then uses Tort method (Tort et al. 2010) to compare, in terms of KL-distance, the ditribution of the amplitude envolope of the high frequencies on the phase of the low frequencies of the signal.    

Inputs: y (your signal, as array), Fs (sampling frequency, integer),  res (resolution of the CWT, integer, default = 140), dec (decimate factor if you want to downsample, integer, default = 10, set as 0 if you do not want to decimate), filt (digital filter, boolean, 1 = yes, 0 = no, default = 1, this will apply a 10 poles lowpass Butterworth IIR forward-backward filter with cutoff at 120Hz), nPerm (integer, default =0, any integer greater than 0 will enable a subroutine that shuffles the amplitudes of the high frequencies on the phase of the recorded signal and computes the MI on the shuffled dataset. This will be used to build a threhosld matrix based on the distribution of MI of the bootstrapped dataset and threosld the observed MI at 1%. Given the computaitonal costs of this operation, this funciton supports parallel computing (see next option)), nJobs (integer, default = 0, any integer greater than 0 enables a subroutine that runs the parallel processing toolbox to spread the MI computation of the shuffled dataset of different cores, nJobs have to ge smaller or equal to the number of the cores of the machine in use. By default it uses the backend “multiprocessor” and “cloupickling” for pickle the variables, those settings can be changed at lines 131-135).  It outputs the Modulation Index, the array of low frequencies (the phase of which is modulating the amplitude of the hgih frequency, default is theta/alpha/beta) the array of the high frequency (the amplitude of which is modulated by the phase of the low frequency array, defulat is gammas), the raw data and Fs (in case decimation and filtering were enabled).


showMI(MI, LfreqE, HfreqE, y, Fs)

Inputs:  modulaiton index, modulating frequency array, modulated frequency array, y (your signal, as array), Fs (sampling frequency, integer). Essentially this routine works with the output of getMI to merge a plot of the raw date and a heatmap of the modulation index. 


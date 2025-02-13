# shift-and-add & grid disentangling, by Tomer Shenar, with contributions from Matthias Fabry & Julia Bodensteiner
# 11.09.2023, V2.0; feel free to contact at T.Shenar@uva.nl or tomer.shenar@gmail.com for questions/inquires
# Algorithm and examples in Gonzales & Levato 2006, A&A, 448, 283; Shenar et al. 2021, A&A, 639, 6; Shenar et al. 2022, A&A, 665, 148
# Function file


import glob
import os
import numpy as np
from astropy.io import ascii
from scipy.interpolate import interp1d
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import scipy.special as sc
from astropy.table import Table
import sys
import pickle
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
legsize = 11
alphaleg = 0.
locleg='lower left'


#Constants
clight = 2.9979E5





                      
## A bunch of functions...

def Normalise(spec, points = []):
    ContIndices = np.array([np.argmin(np.abs(spec[:,0] - points[i])) for i in np.arange(len(points))])  
    Contfluxes = np.array([np.average(spec[:,1][max(ContIndices[i] - 7, 0): min(ContIndices[i] + 7, len(spec[:,0]))]) for i in np.arange(len(ContIndices))])
    #print spec[:,0][ContIndices]
    #print Contfluxes
    #dasd
    ContSpline = interp1d(spec[:,0][ContIndices], Contfluxes, bounds_error = False, fill_value = 'extrapolate')(spec[:,0])
    return np.array([spec[:,0], spec[:,1] / ContSpline]).T 


    print(("Data written to %s" % outfilename))


def Cosclean(Spec, thold=6, cosize=10, ForbiddenRanges = [[3970, 3975], [4020, 4030], [4103., 4108.], [4342, 4347], [4365, 4369],  [4391, 4393], [4470, 4477.]]):
    for itr in range(10):
        waves = np.copy(Spec[:,0])
        fluxes = np.copy(Spec[:,1])
        for wrange in ForbiddenRanges:
            if 'WaveCond' not in locals():
                WaveCond = (waves > wrange[0]) * (waves < wrange[1])
            else:
                WaveCond += (waves > wrange[0]) * (waves < wrange[1])
        fluxes[WaveCond] = 1.        
        fluxdiff = np.append(0, np.diff(fluxes))
        sigma =  thold*np.mean(np.absolute(fluxdiff))
#Find points whose gradients are larger than thold*average
        gradient_condition = (np.absolute(fluxdiff) > sigma)     
#Weak point are of no interest
        flux_condition = fluxes > 1.0
        posgrad = fluxdiff>0
        flagged_args = np.where(gradient_condition & flux_condition)[0]
        if (not len(flagged_args)):
            print("There are no cosmics detected with the given threshhold and size")
            return Spec
        blimit = 0
        N = len(waves)
        for i in flagged_args:
            if waves[i] < blimit or i<cosize or i>N-cosize:
                continue
            cosmic = fluxes[i-cosize:i+cosize+1]    
            if posgrad[i]:
                if np.any(cosmic[:cosize] - fluxes[i] > 0):
                    continue
            else:
                if np.any(cosmic[cosize+1:] - fluxes[i] > 0):
                    continue    
            ipeak = i - cosize + np.argmax(fluxes[i-cosize:i+cosize+1])
            fpeak = fluxes[ipeak]
            cosmic = fluxes[ipeak - cosize : ipeak + cosize + 1]
            cosb = cosmic[:cosize + 1]
            cosr = cosmic[cosize:]
            fmeadb = np.mean(cosb)
            fmeadr = np.mean(cosr)
            cosbdiff = np.append(np.diff(cosb), 0)
            cosrdiff = np.append(0, np.diff(cosr))
            sigmab = np.mean(np.absolute(cosbdiff))
            sigmar = np.mean(np.absolute(cosrdiff))       
            condsmallb = cosb - fmeadb < 0.1*(fpeak - fmeadb)
            condsmallr = cosr - fmeadr < 0.1*(fpeak - fmeadr)
            argb = np.where((np.roll(cosbdiff,-1) > sigmab) & (condsmallb) & (cosb > 0.5))[0]
            argr = np.where((cosrdiff < -sigmar) & (condsmallr) & (cosr > 0.5))[0]
            if len(argb) == 0  or len(argr) == 0:
                continue
            argb = ipeak - cosize + argb[-1]
            argr = ipeak + argr[0]
            if abs(fluxes[argb] - fpeak) < sigmab or abs(fluxes[argr] - fpeak) < sigmar:
                continue
            Spec[argb:argr+1,1] = np.interp(waves[argb:argr+1], [waves[argb], waves[argr]], [fluxes[argb], fluxes[argr]]) 
            blimit = waves[argr]
    return Spec


    
# Solves the Kepler equation
def Kepler(E, M, ecc):
      E2 = (M - ecc*(E*np.cos(E) - np.sin(E))) / (1. - ecc*np.cos(E))
      eps = np.abs(E2 - E) 
      if np.all(eps < 1E-10):
            return E2
      else:
            return Kepler(E2, M, ecc)




# Given true anomaly nu and parameters, 
def v1andv2(nu, Orbital_Params, Config = 'in'):
      P, T0, ecc, omega, Gamma, K1, K2, Pout, T0Out, eccOut, omegaOut, KOut, P2, T0_2, e2, omega2, K3, K4  = list(Orbital_Params.values()) 
      #print(P, Pout, ecc, eccOut, K1, K2, KOut, omega, omegaOut)
      Omegapi = omega/180. * np.pi  
      OmegapiOut = omegaOut/180. * np.pi 
      v1 = Gamma + K1*(np.cos(Omegapi + nu) + ecc* np.cos(Omegapi))   
      v2 = Gamma - K2*(np.cos(Omegapi + nu) + ecc* np.cos(Omegapi))    
      if Config=='in':
        return v1, v2
      elif Config=='Tertiary':
        K_bin = ( P/Pout * ( (1 - ecc**2)/(1 - eccOut**2) )**(3./2.) * (K1 + K2)**3 / KOut )**0.5 - KOut
        if K_bin < 0:
          print("WARNING: this configuration is not physical for inc_out = inc_in!!!")
        vOut = Gamma + KOut*(np.cos(OmegapiOut + nu) + eccOut* np.cos(OmegapiOut))   
        vBin = Gamma - K_bin*(np.cos(OmegapiOut + nu) + eccOut* np.cos(OmegapiOut))     
        return vOut, vBin
    
    



#Calculate K corresponding to chi2 minimum by fitting parabola to minimum region + confidence interval (default: 1sig=68%)    
# redchi2: array of reduced chi2, i.e chi2/DoF (DoF = Degrees of Freedom)
# nu = DoF
def Chi2con(redchi2, nu, Kscomp, Rangestr, P1=0.68, comp='secondary', ParbSize = 3):    
#First fit parabola to chi2 distribution to find minimum (indmin-ParbSize, indmin+ParbSize, default ParbSize=3)
    indmin = np.argmin(redchi2)    
    i1 = max(0, indmin-ParbSize)
    i2 = min(len(Kscomp)-1, indmin+ParbSize)
    a,b,c = np.polyfit(Kscomp[i1:i2], redchi2[i1:i2], 2)      
    if a < 0:
        print("Could not fit sensible parabola (a < 0); try changing the fitting range by increasing/decreasing ParbSize argument")
        plt.scatter(Kscomp, redchi2)
        plt.show()
        return 0., Kscomp[indmin], 0.
    ParbMin = c - b**2/4./a
# Now compute non-reduced, normalised chi2 distribution  
    chi2 = redchi2 * nu /  ParbMin
# The probability distribution of (non-reduced) chi2 peaks at nu. Compute array around this value    
    xarr = np.arange(nu/10., nu*10, 1)
# The cumulative probability distribution of chi^2 (i.e., Prob(chi2)<x) with nu DoFs is the regularised incomplete Gamma function Gamma(nu/2, x/2)
# Hence, we look for the index and x value at which Prob(chi2) < P1
    ys1 = sc.gammainc(nu/2., xarr/2) - P1
    minarg1 = np.argmin(np.abs(ys1))
# This is the chi2 value corresponding to P1 (typically 1-sigma)    
    chi2P1 = xarr[minarg1]/nu   
# Now fit parabola to reduced chi2, after normalisation:
    a,b,c = np.polyfit(Kscomp[i1:i2], chi2[i1:i2]/nu, 2)  
    chi2fine = np.arange(Kscomp[i1], Kscomp[i2], 0.01)
    parb = a*chi2fine**2 + b*chi2fine  + c
    plt.scatter(Kscomp, chi2/nu*ParbMin)
    plt.plot(chi2fine, parb*ParbMin, color='orange')    
    K2min = -b/2./a
    K2err = K2min - (-b - np.sqrt(b**2 - 4*a*(c-chi2P1))) / 2./a
    plt.plot([Kscomp[0], Kscomp[-1]], [chi2P1*ParbMin, chi2P1*ParbMin], color='red', label=r'1$\sigma$ contour')    
    plt.ylabel(r'Reduced $\chi^2$')        
    plt.legend()
    if comp=='secondary':
        plt.xlabel(r'$K_2$ [km/s]')
        np.savetxt('disentangled/' + Rangestr + '_' + 'grid_dis_K2.txt', np.c_[Kscomp, redchi2], header='#1sigma = ' + str(chi2P1*ParbMin))     
        plt.savefig('disentangled/' +  Rangestr +  '_Grid_disentangling_K2.pdf', bbox_inches='tight')        
    elif comp=='primary':
        plt.xlabel(r'$K_1$ [km/s]')        
        np.savetxt('disentangled/' + Rangestr + '_' + 'grid_dis_K1.txt', np.c_[Kscomp, redchi2], header='#1sigma = ' + str(chi2P1*ParbMin))    
        plt.savefig('disentangled/' + Rangestr +  '_Grid_disentangling_K1.pdf', bbox_inches='tight')        
    elif comp=='tertiary':
        plt.xlabel(r'$K_3$ [km/s]')        
        np.savetxt('disentangled/' + Rangestr + '_' + 'grid_dis_K3.txt', np.c_[Kscomp, redchi2], header='#1sigma = ' + str(chi2P1*ParbMin))    
        plt.savefig('disentangled/' + Rangestr +  '_Grid_disentangling_K3.pdf', bbox_inches='tight')       
    else:
        sys.exit("type of companion not defined -- error in subroutine Chi2con...")        
    plt.show()    
    return chi2P1, K2min, K2err

# Reead HERMES data
def read_HERMES(infile):
    print(("%s: Input file is a HERMES file." % infile))
    header = fits.getheader(infile)
    # for files with standard wavelegth array
    if ((header['CTYPE1'] == 'WAVELENGTH') or (header['CTYPE1'] == 'AWAV')):
        flux = fits.getdata(infile)
        crval = header['CRVAL1']
        cdelt = header['CDELT1']
        naxis1 = header['NAXIS1']
        wave = crval + np.arange(0, naxis1) * cdelt

    # for files that are given in logarithmic wl array
    if (header['CTYPE1'] == 'log(wavelength)'):
        flux = fits.getdata(infile)
        crval = header['CRVAL1']
        cdelt = header['CDELT1']
        naxis1 = header['NAXIS1']
        wave = np.exp(crval + np.arange(0, naxis1)*cdelt)

    else:
        print("Could not read in HERMES fits file - unknown file type.")
        sys.exit()
    flux = np.nan_to_num(flux, 1.)
    return np.array([wave, flux]).T

#Ensures that arr values where arr > lim  are set to 0 in domains specified in Poslim array
def Limit(waves, arr, lim, Poslim):
    for i, Range in enumerate(Poslim):
        if 'PosCond' not in locals():
            PosCond = (Poslim[i][0] < waves)  * (waves < Poslim[i][1])
        else:
            PosCond += (Poslim[i][0] < waves)  * (waves < Poslim[i][1])
    NegCond = np.logical_not(PosCond)    
    arr[(arr > lim) * NegCond] = 0.
    return arr




# Shrinks the wavelength domain on which obs-mod is calculated to avoid edge issues
def Reduce_Waves(waves, nusdata, Orbital_Params, K1max, K2max):
    Orbital_Params_max = Orbital_Params.copy()
    Orbital_Params_max['K1'] = K1max
    Orbital_Params_max['K2'] = K2max    
    vrA, vrB = v1andv2(nusdata, Orbital_Params_max)
    Inds = np.where(np.diff(waves) > waves[1]-waves[0])
    WaveCalcCond = waves < 0
    if len(Inds) == 0:
        LamMin = waves[0]*(1. +  max(max(vrA), max(vrB))/clight)
        LamMax = waves[-1]*(1. +  min(min(vrA), min(vrB))/clight)
        WaveCalcCond = (waves > LamMin) * (waves < LamMax)
    else:
        Inds = np.append(0, Inds[0])
        Inds = np.append(Inds, len(waves)-1)
        for j in np.arange(len(Inds)-1):
            LamMin = waves[Inds[j]+1]*(1. +  max(max(vrA), max(vrB))/clight)
            LamMax = waves[Inds[j+1]]*(1. +  min(min(vrA), min(vrB))/clight)
            WaveCalcCond = WaveCalcCond + (waves > LamMin) * (waves < LamMax)
##Reduce nebular lines regions:
    #if NebLineCutchi2 == True:
        #for wrange in PosLimCondC:
            #if 'NebCond' in locals():
                #NebCond += (waves > wrange[0]*(1+Gamma/clight) ) * (waves < wrange[1]*(1+Gamma/clight)  )
            #else:
                #NebCond = (waves > wrange[0]*(1+Gamma/clight)) * (waves < wrange[1]*(1+Gamma/clight))
        #WaveCalcCond *= ~NebCond
    return WaveCalcCond


def Prepare_Plot_Extremes(fig, axes, waves, SpecShiftArr, Nebshift, ObsSpec, specsum, specname, PhiArr, MJD, pltExtyMin, pltExtyMax, StarName, Rangestr, KNowArr, PanelCount, linewidth=3, NebLines=False, Panel=0, linewidExt=3, ExtremesFigSize=(7,8), ResizeExtremesFac_min = 1., ResizeExtremesFac_max = 1., TitleFileName=False, DisLine=None, ObsCol = 'blue', PrimCol='black', SecCol='green', TerCol='purple', SumCol='red', NebCol='pink' ):
    axes[Panel].plot(waves, SpecShiftArr[0], label='Prim Dis.', color=PrimCol, linestyle = 'dotted', linewidth=linewidExt)
    axes[Panel].plot(waves, SpecShiftArr[1], label='Sec. Dis.', color=SecCol, linewidth=linewidExt)   
    if NebLines:
        axes[Panel].plot(waves, Nebshift, label='Nebular Dis.', color=NebCol)               
    if len(SpecShiftArr)==3:
        axes[Panel].plot(waves, SpecShiftArr[2], label='Ter. Dis.', color=TerCol, linewidth=linewidExt)           
    if len(SpecShiftArr)==2:
        axes[Panel].plot(waves, ObsSpec, color=ObsCol, label=str(round(MJD,1)) + r', $\varphi=$' + str(round(PhiArr[0], 2)))
    elif len(SpecShiftArr)==3:
        axes[Panel].plot(waves, ObsSpec, color=ObsCol, label=str(round(MJD,1)) + r', $\varphi_{\rm in}=$' + str(round(PhiArr[0], 2)) + r', $\varphi_{\rm Out}=$' + str(round(PhiArr[1], 2))    )        
    axes[Panel].plot(waves, specsum, label='Sum Dis.', color=SumCol, linestyle = '--', linewidth=linewidExt)    
    if TitleFileName == True:
        axes[Panel].set_title(specname.split('/')[-1])
    axes[Panel].legend(prop={'size': legsize}, loc=locleg, framealpha=alphaleg)
    axes[Panel].set_ylabel('Normalised flux')
    axes[Panel].set_ylim(ResizeExtremesFac_min*pltExtyMin, ResizeExtremesFac_max*pltExtyMax)
    DiffMajor = max(int((waves[-1] - waves[0])/3), 1)
    axes[Panel].xaxis.set_major_locator(MultipleLocator(DiffMajor))
    axes[Panel].xaxis.set_minor_locator(MultipleLocator(DiffMajor/5.))
    axes[Panel].yaxis.set_minor_locator(MultipleLocator(.01))   
    axes[Panel].plot([waves[0], waves[-1]], [0., 0.], color='black')   
    #plt.show()
    if DisLine != None and Panel==0:
        axes[Panel].text(0.75, 0.1, DisLine, size=20, horizontalalignment='center',  verticalalignment='center', transform = axes[Panel].transAxes)    
    PanelCount +=1
    if len(SpecShiftArr)==2:
        FinPanelNum = 2
    elif len(SpecShiftArr)==3:
        FinPanelNum = 4
    if PanelCount==FinPanelNum:
        axes[Panel].set_xlabel(r'Wavelength$\,[\AA]$')
        plt.tight_layout()    
        if len(SpecShiftArr)==2:
            try:
                NameFile =  'disentangled/' + StarName + '_' + Rangestr + '_Extremes_' + str(np.round(KNowArr[0])) + '_' + str(np.round(KNowArr[1])) + '.pdf'   
                plt.savefig(NameFile, bbox_inches='tight')                                       
            except:
                NameFile =  'disentangled/' + StarName + '_' + Rangestr + '_Extremes_' + str(np.round(Orbital_Params['K1'])) + '_' + str(np.round(K2s[kcount])) + '.pdf'   
                plt.savefig(NameFile, bbox_inches='tight')   
        elif len(SpecShiftArr)==3:
            try:
                NameFile =  'disentangled/' + StarName + '_' + Rangestr + '_Extremes_' + str(np.round(KNowArr[0])) + '_' + str(np.round(KNowArr[1])) + '_' + str(np.round(KNowArr[2])) +  '.pdf'   
                plt.savefig(NameFile, bbox_inches='tight')                                       
            except:
                sys.error("No idea why this should fail... Maybe you do? Existing 'Prepare_Plot_Extremes' routine")           
        plt.show()    
    return PanelCount


# Calculate difference 
def CalcDiffs(DisSpecVector, vrA, vrB,  waves, ObsSpecs, nusdata,  Orbital_Params, K1s, K2s, MJDs, phis, specnames, Rangestr, StarName, ScalingNeb, Resid=False, Reduce=False, ShowItr=False, PLOTEXTREMES=False, PLOTFITS=False, kcount_extremeplot=0,  linewidExt=3, NebLines = False, NebFac=1, S2NpixelRange=5, kcount_usr=0, ExtremesFigSize=(7,8), ResizeExtremesFac_max=1, ResizeExtremesFac_min=1, TitleFileName=False, DisLine=None):
    global kcount, k1, k2, DoFs, K1now, K2now
    RVExtmaxInd, RVExtminInd = np.argmax(vrA), np.argmin(vrA)     
    A, B, NebSpec = DisSpecVector
    if Resid:
        Residuals = []    
    if 'DoFs' not in globals():
        DoFs=1
    WaveCalcCond = Reduce_Waves(waves, nusdata, Orbital_Params, K1s[-1], K2s[-1])
    Sum = 0     
    if 'kcount' not in globals():
        kcount = 0
    if PLOTEXTREMES: 
        plotminyarr=[]
        plotmaxyarr=[]
        for ind in np.arange(len(ObsSpecs)):
            plotminyarr.append(np.amin(interp1d(ObsSpecs[ind][:,0], ObsSpecs[ind][:,1]-1,bounds_error=False, fill_value=0.)(waves[WaveCalcCond])))
            plotmaxyarr.append(np.amax(interp1d(ObsSpecs[ind][:,0], ObsSpecs[ind][:,1]-1,bounds_error=False, fill_value=0.)(waves[WaveCalcCond])))
        plotminyarr.append(min(A))
        plotminyarr.append(min(B)) 
        plotmaxyarr.append(max(A))
        plotmaxyarr.append(max(B))        
        pltExtyMin = min(plotminyarr)*1.1
        pltExtyMax = max(plotmaxyarr)*1.1
    ExtremePlotDefined = False           
    for ind in np.arange(len(ObsSpecs)):
        vA = vrA[ind]/clight
        vB = vrB[ind]/clight
        Facshift1 = np.sqrt( (1 + vA) / (1 - vA))
        Facshift2 = np.sqrt( (1 + vB) / (1 - vB))
        Ashift = interp1d(waves*Facshift1, A,bounds_error=False, fill_value=0.)(waves[WaveCalcCond])   
        Bshift = interp1d(waves*Facshift2, B,bounds_error=False, fill_value=0.)(waves[WaveCalcCond]) 
        ObsSpec = interp1d(ObsSpecs[ind][:,0], ObsSpecs[ind][:,1]-1,bounds_error=False, fill_value=0.)(waves[WaveCalcCond])   
        if NebLines:
            vNeb = 0.                    
            Nebshift = NebFac*ScalingNeb[ind]*interp1d(waves, NebSpec,bounds_error=False, fill_value=0.)(waves[WaveCalcCond])   
            specsum = Ashift + Bshift + Nebshift   
        else:      
            Nebshift = np.zeros(len(Ashift))
            specsum = Ashift + Bshift              
        sigma =  (np.std(ObsSpec[:S2NpixelRange]) +  np.std(ObsSpec[-S2NpixelRange:]))/2.

        if Resid:
            Residuals.append(ObsSpec - specsum)        
        Sum +=np.sum( (ObsSpec - specsum)**2/sigma**2)   
        if PLOTFITS:
# User needs to change kcount==1 condition if a specific K2 is desired for plotting.             
            if kcount==kcount_usr:      
                plt.plot(waves[WaveCalcCond], specsum, label='sum')
                plt.plot(waves[WaveCalcCond], Ashift, label='A')
                plt.plot(waves[WaveCalcCond], Bshift, label='B')   
                if NebLines:
                    plt.plot(waves[WaveCalcCond], Nebshift, label='Neb')   
                plt.plot(waves[WaveCalcCond], ObsSpec, label=specnames[ind].split('/')[-1]  + r', $\varphi=$' + str(round(phis[ind], 2)))
                plt.plot([waves[WaveCalcCond][0], waves[WaveCalcCond][-1]], [0., 0.], color='black')
                plt.legend()
                plt.show()   
        if PLOTEXTREMES:
            if kcount==kcount_extremeplot:
                if not ExtremePlotDefined:
                    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=ExtremesFigSize)  
                    ExtremePlotDefined = True      
                    PanelCount = 0                 
                ShiftSpecArr = [Ashift, Bshift]
                PhiArr = [phis[ind]]
                KNowArr = [K1now, K2now]
                if ind==min(RVExtminInd, RVExtmaxInd):
                    PanelCount = Prepare_Plot_Extremes(fig, axes, waves[WaveCalcCond], ShiftSpecArr, Nebshift, ObsSpec, specsum, specnames[ind], PhiArr, MJDs[ind], pltExtyMin, pltExtyMax, StarName, Rangestr, KNowArr, PanelCount, linewidth=linewidExt, NebLines=NebLines, ExtremesFigSize=ExtremesFigSize, Panel=0, ResizeExtremesFac_max = ResizeExtremesFac_max, ResizeExtremesFac_min = ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)                   
                elif ind==max(RVExtminInd, RVExtmaxInd):
                    PanelCount = Prepare_Plot_Extremes(fig, axes, waves[WaveCalcCond], ShiftSpecArr, Nebshift, ObsSpec, specsum, specnames[ind], PhiArr, MJDs[ind], pltExtyMin, pltExtyMax, StarName, Rangestr, KNowArr, PanelCount, linewidth=linewidExt, NebLines=NebLines,  ExtremesFigSize=ExtremesFigSize, Panel=1, ResizeExtremesFac_max = ResizeExtremesFac_max, ResizeExtremesFac_min = ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)  
    print("kcount:", kcount)   
    if DoFs==0:
        try:
            DoFs += len(waves[WaveCalcCond]) * (len(ObsSpecs)-2)
        except:
            pass
    kcount+=1
    if not ShowItr:
        print("chi2:",  Sum/ (len(waves[WaveCalcCond])* len(ObsSpecs) - 2))
    if Resid:
        return Residuals
    if Reduce:
        return Sum/ (len(waves[WaveCalcCond]) * len(ObsSpecs) - 1)
    else:
        return Sum
    
    

# Calculate difference 
def CalcDiffs3D(DisSpecVector, vrA, vrB, vrOut, vrBin, waves, ObsSpecs, nusdata, nusdataOut,  Orbital_Params, K1s, K2s, KOuts, MJDs, phis, phisOut, specnames, Rangestr, StarName, ScalingNeb, Resid=False, Reduce=False, ShowItr=False, PLOTEXTREMES=False, PLOTFITS=False, kcount_extremeplot=0,  linewidExt=3, NebLines = False, NebFac=1, S2NpixelRange=5, kcount_usr=0, ExtremesFigSize=(7,8), ResizeExtremesFac_max=1, ResizeExtremesFac_min=1, TitleFileName=False, DisLine=None):
    global kcount, k1, k2, kOut, DoFs, K1now, K2now, K3now
    RVExtmaxInd, RVExtminInd = np.argmax(vrA), np.argmin(vrA)     
    RVExtmaxInd_Out, RVExtminInd_Out = np.argmax(vrOut), np.argmin(vrOut)     
    A, B, C, NebSpec = DisSpecVector
    if Resid:
        Residuals = []    
    if 'DoFs' not in globals():
        DoFs=1
    WaveCalcCond = Reduce_Waves(waves, nusdata, Orbital_Params, K1s[-1], K2s[-1])
    Sum = 0     
    if 'kcount' not in globals():
        kcount = 0
    if PLOTEXTREMES: 
        plotminyarr=[]
        plotmaxyarr=[]
        for ind in np.arange(len(ObsSpecs)):
            plotminyarr.append(np.amin(interp1d(ObsSpecs[ind][:,0], ObsSpecs[ind][:,1]-1,bounds_error=False, fill_value=0.)(waves[WaveCalcCond])))
            plotmaxyarr.append(np.amax(interp1d(ObsSpecs[ind][:,0], ObsSpecs[ind][:,1]-1,bounds_error=False, fill_value=0.)(waves[WaveCalcCond])))
        plotminyarr.append(min(A))
        plotminyarr.append(min(B)) 
        plotmaxyarr.append(max(A))
        plotmaxyarr.append(max(B))        
        pltExtyMin = min(plotminyarr)*1.1
        pltExtyMax = max(plotmaxyarr)*1.1
    ExtremePlotDefined = False        
    for ind in np.arange(len(ObsSpecs)):
        vA = (vrA[ind] + vrBin[ind])/clight
        vB = (vrB[ind] + vrBin[ind])/clight
        vOut = vrOut[ind]/clight
        Facshift1 = np.sqrt( (1 + vA) / (1 - vA))
        Facshift2 = np.sqrt( (1 + vB) / (1 - vB))
        FacshiftOut = np.sqrt( (1 + vOut) / (1 - vOut))
        Ashift = interp1d(waves*Facshift1, A,bounds_error=False, fill_value=0.)(waves[WaveCalcCond])   
        Bshift = interp1d(waves*Facshift2, B,bounds_error=False, fill_value=0.)(waves[WaveCalcCond]) 
        Cshift = interp1d(waves*FacshiftOut, C,bounds_error=False, fill_value=0.)(waves[WaveCalcCond]) 
        ObsSpec = interp1d(ObsSpecs[ind][:,0], ObsSpecs[ind][:,1]-1,bounds_error=False, fill_value=0.)(waves[WaveCalcCond])   
        if NebLines:
            vNeb = 0.                    
            Nebshift = NebFac*ScalingNeb[ind]*interp1d(waves, NebSpec,bounds_error=False, fill_value=0.)(waves[WaveCalcCond])   
            specsum = Ashift + Bshift + Cshift +  Nebshift   
        else:      
            Nebshift = np.zeros(len(Ashift))
        specsum = Ashift + Bshift + Cshift +  Nebshift             
        sigma =  (np.std(ObsSpec[:S2NpixelRange]) +  np.std(ObsSpec[-S2NpixelRange:]))/2.
        if Resid:
            Residuals.append(ObsSpec - specsum)        
        Sum +=np.sum( (ObsSpec - specsum)**2/sigma**2)   
        if PLOTFITS:
# User needs to change kcount==1 condition if a specific K2 is desired for plotting.             
            if kcount==kcount_usr:      
                plt.plot(waves[WaveCalcCond], specsum, label='sum')
                plt.plot(waves[WaveCalcCond], Ashift, label='A')
                plt.plot(waves[WaveCalcCond], Bshift, label='B')   
                plt.plot(waves[WaveCalcCond], Cshift, label='C')   
                if NebLines:
                    plt.plot(waves[WaveCalcCond], Nebshift, label='Neb')   
                plt.plot(waves[WaveCalcCond], ObsSpec, label=specnames[ind].split('/')[-1]  + r', $\varphi=$' + str(round(phis[ind], 2)))
                plt.plot([waves[WaveCalcCond][0], waves[WaveCalcCond][-1]], [0., 0.], color='black')
                plt.legend()
                plt.show()   
        if PLOTEXTREMES:
            if kcount==kcount_extremeplot:
                if not ExtremePlotDefined:
                    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=ExtremesFigSize)  
                    ExtremePlotDefined = True
                    PanelCount = 0
                ShiftSpecArr = [Ashift, Bshift, Cshift]
                PhiArr = [phis[ind], phisOut[ind]]
                KNowArr = [K1now, K2now, K3now]    
                if ind==RVExtminInd:
                    PanelCount = Prepare_Plot_Extremes(fig, axes, waves[WaveCalcCond], ShiftSpecArr, Nebshift, ObsSpec, specsum, specnames[ind], PhiArr,  MJDs[ind], pltExtyMin, pltExtyMax, StarName, Rangestr, KNowArr, PanelCount, linewidth=linewidExt, NebLines=NebLines, ExtremesFigSize=ExtremesFigSize, Panel=(0,0), ResizeExtremesFac_max = ResizeExtremesFac_max, ResizeExtremesFac_min = ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)                   
                if ind==RVExtmaxInd:
                    PanelCount = Prepare_Plot_Extremes(fig, axes, waves[WaveCalcCond], ShiftSpecArr, Nebshift, ObsSpec, specsum, specnames[ind], PhiArr, MJDs[ind], pltExtyMin, pltExtyMax, StarName, Rangestr, KNowArr, PanelCount, linewidth=linewidExt, NebLines=NebLines,  ExtremesFigSize=ExtremesFigSize, Panel=(1,0), ResizeExtremesFac_max = ResizeExtremesFac_max, ResizeExtremesFac_min = ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)  
                if ind==RVExtmaxInd_Out:
                    PanelCount = Prepare_Plot_Extremes(fig, axes, waves[WaveCalcCond], ShiftSpecArr, Nebshift, ObsSpec, specsum, specnames[ind], PhiArr,  MJDs[ind], pltExtyMin, pltExtyMax, StarName, Rangestr, KNowArr, PanelCount, linewidth=linewidExt, NebLines=NebLines,  ExtremesFigSize=ExtremesFigSize, Panel=(0,1), ResizeExtremesFac_max = ResizeExtremesFac_max, ResizeExtremesFac_min = ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)  
                if ind==RVExtminInd_Out:
                    PanelCount = Prepare_Plot_Extremes(fig, axes, waves[WaveCalcCond], ShiftSpecArr, Nebshift, ObsSpec, specsum, specnames[ind], PhiArr, MJDs[ind], pltExtyMin, pltExtyMax, StarName, Rangestr,KNowArr, PanelCount, linewidth=linewidExt, NebLines=NebLines,  ExtremesFigSize=ExtremesFigSize, Panel=(1,1), ResizeExtremesFac_max = ResizeExtremesFac_max, ResizeExtremesFac_min = ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)  
    print("kcount:", kcount)   
    if DoFs==0:
        try:
            DoFs += len(waves[WaveCalcCond]) * (len(ObsSpecs)-3)
        except:
            pass
    kcount+=1
    if not ShowItr:
        print("chi2:",  Sum/ (len(waves[WaveCalcCond])* len(ObsSpecs) - 3))
    if Resid:
        return Residuals
    if Reduce:
        return Sum/ (len(waves[WaveCalcCond]) * len(ObsSpecs) - 1)
    else:
        return Sum    



# Documentation: 
# "B" = array, initial guess for flux of "secondary"
# vrads1, vrads2 = RVs of primary, secondary
# waves: wavelength grid on which disentanglement should take place
#NOTE: If initial guess for primary is preferred, roles of primary should change, i.e., one should call:
# disentangle(Aini, vrads2, vrads1, waves)
# Resid --> returns array of residual spectra between obs and dis1+dis2
# Reduce --> Returns reduced chi2
def disentangle(B, vrads1, vrads2, waves, ObsSpecs, weights, StrictNeg, PosLimCond, Poslimall,  nusdata, Orbital_Params, K1s, K2s, MJDs, phis, specnames,  Rangestr, StarName, ScalingNeb, NebSpec, Resid=False, Reduce=False, ShowItr=False, Once=False,  InterKind='linear', itrnumlim=100, PLOTCONV=False, PLOTITR=False,   PLOTEXTREMES=False, PLOTFITS=False, kcount_extremeplot=0, linewidExt=3,  N_Iteration_Plot=50, NebLines = False, NebFac=1, kcount_usr=0, ExtremesFigSize=(7,8), ResizeExtremesFac_max=1., ResizeExtremesFac_min=1, TitleFileName=False, DisLine=None):    
    global kcount, k1, k2, DoFs, kcount, K1now, K2now
    StrictNegA, StrictNegB, StrictNegC, StrictNegD = StrictNeg
# If convergence plot: allow spectra to be positive for sensible convergence plot:
    if PLOTCONV:
        StrictNegA, StrictNegB, StrictNegC, StrictNegD = False, False, False, False    
    PosLimCondA, PosLimCondB, PosLimCondC, PosLimCondD, PosLimCondNeb = PosLimCond
    if (not Once):
        try:
            K1now = K1s[k1]
            K2now = K2s[k2]   
            print("Disentangeling..... K1, K2=", K1now, K2now)
        except:
            pass
    else:
        k1, k2 = 0, 0
    Facshift1 = np.sqrt( (1 + vrads1/clight) / (1 - vrads1/clight))
    Facshift2 = np.sqrt( (1 + vrads2/clight) / (1 - vrads2/clight))
    Ss1 = [interp1d(ObsSpecs[i][:, 0] /Facshift1[i], ObsSpecs[i][:, 1]-1.,
                      bounds_error=False, fill_value=0., kind=InterKind)(waves) for i in np.arange(len(ObsSpecs))] 
    Ss2 = [interp1d(ObsSpecs[i][:, 0] /Facshift2[i], ObsSpecs[i][:, 1]-1.,
                      bounds_error=False, fill_value=0., kind=InterKind)(waves) for i in np.arange(len(ObsSpecs))]      
#Frame of Refernce star 1:    
    WavesBA = np.array([waves * Facshift2[i] /Facshift1[i] for i in np.arange(len(vrads1))])  
#Frame of Refernce star 2:        
    WavesAB = np.array([waves * Facshift1[i] /Facshift2[i] for i in np.arange(len(vrads1))])  
    if NebLines:
        FacshiftNeb = np.ones(len(vrads1))
        SsNeb = [interp1d(ObsSpecs[i][:, 0] / FacshiftNeb[i], ObsSpecs[i][:, 1]-1.,
                      bounds_error=False, fill_value=0., kind=InterKind)(waves) for i in np.arange(len(ObsSpecs))]             
        WavesNebA = np.array([waves *FacshiftNeb[i]/Facshift1[i] for i in np.arange(len(vrads1))])  
        WavesNebB = np.array([waves * FacshiftNeb[i] /Facshift2[i] for i in np.arange(len(vrads1))])         
        WavesANeb = np.array([waves *Facshift1[i] /FacshiftNeb[i] for i in np.arange(len(vrads1))])  
        WavesBNeb = np.array([waves * Facshift2[i] /FacshiftNeb[i] for i in np.arange(len(vrads1))])   
    itr = 0
    while itr<itrnumlim:
        itr+=1
        BAshifts = [interp1d(WavesBA[i], B,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]    
        if NebLines:
            NebAshifts = [interp1d(WavesNebA[i], NebSpec,bounds_error=False, fill_value=0., kind=InterKind) for i in np.arange(len(ObsSpecs))]    
            SpecMean = np.sum(np.array([[i]*(Ss1[i] - BAshifts[i](waves) - NebFac*ScalingNeb[i]*NebAshifts[i](waves)) for i in np.arange(len(Ss1))]), axis=0)      
        else:
            SpecMean = np.sum(np.array([weights[i]*(Ss1[i] - BAshifts[i](waves) ) for i in np.arange(len(Ss1))]), axis=0)          
        Anew = interp1d(waves, SpecMean, bounds_error=False, fill_value=0., kind=InterKind)(waves)  
        if StrictNegA:
            Anew = Limit(waves, Anew, Poslimall, PosLimCondA)        
        if 'A' in locals():
            Epsnew = np.amax((A - Anew)**2)        
        else:
            Epsnew = 0.
        A = np.copy(Anew)   
        ABshifts = [interp1d(WavesAB[i], A,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]     
        if NebLines:
            NebBshifts = [interp1d(WavesNebB[i], NebSpec,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]     
            SpecMean = np.sum(np.array([weights[i]*(Ss2[i] - ABshifts[i](waves) - NebFac*ScalingNeb[i]*NebBshifts[i](waves)) for i in np.arange(len(Ss1))]), axis=0)
        else:
            SpecMean = np.sum(np.array([weights[i]*(Ss2[i] - ABshifts[i](waves)) for i in np.arange(len(Ss1))]), axis=0)                    
        Bnew = interp1d(waves, SpecMean, bounds_error=False, fill_value=0., kind=InterKind)(waves) 
        if StrictNegB:
            Bnew = Limit(waves, Bnew, Poslimall, PosLimCondB)            
        Epsnew = max(Epsnew, np.sum((B - Bnew)**2))
        B = Bnew 
        if NebLines:
            ANebshifts = [interp1d(WavesANeb[i], A,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]     
            BNebshifts = [interp1d(WavesBNeb[i], B,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]                 
            SpecMean = np.sum(np.array([weights[i]*Limit(waves, SsNeb[i] - ANebshifts[i](waves) - BNebshifts[i](waves), Poslimall, PosLimCondNeb) for i in np.arange(len(Ss1))]), axis=0)   
            NebSpecnew = interp1d(waves, SpecMean, bounds_error=False, fill_value=0., kind=InterKind)(waves)
            NebSpecnew[NebSpecnew<Poslimall[0]] = 0.
            Epsnew = max(Epsnew, np.sum(np.abs(NebSpec - NebSpecnew))) 
            NebSpec = NebSpecnew           
        if ShowItr:
            if (itr/itrnumlim*100)%10==0:
                print("Finished " + str(itr) + " out of " + str(itrnumlim) + " iterations (" +str(np.round(itr/itrnumlim*100.,3)) + "%)")
                #print("Convergence Epsilon:", Epsnew)
        if PLOTCONV:
            plt.scatter(itr, np.log10(Epsnew), color='blue')
        if PLOTITR:
            if itr%N_Iteration_Plot==0:
                plt.plot(waves, A, label=itr)
                plt.plot(waves, B, label=itr)
                if NebLines:
                    plt.plot(waves, NebSpec, label=itr)
    print("Finished after ", itr, " iterations")
    if PLOTCONV or PLOTITR:
        if PLOTITR:
            plt.ylabel('Normalised flux')
            plt.xlabel('Wavelength')
            plt.legend()
        if PLOTCONV:
            plt.ylabel('log(Eps)')
            plt.xlabel('iteration number')
        plt.show()  
    DisSpecVector = np.array([A, B , NebSpec])
    return DisSpecVector+1.,  CalcDiffs(DisSpecVector, vrads1,  vrads2,  waves, ObsSpecs, nusdata, Orbital_Params, K1s, K2s, MJDs, phis, specnames, Rangestr, StarName, ScalingNeb, Resid=Resid, Reduce=Reduce, ShowItr=ShowItr,   PLOTEXTREMES=PLOTEXTREMES, PLOTFITS=PLOTFITS, kcount_extremeplot=kcount_extremeplot,  linewidExt=linewidExt, NebLines = NebLines, NebFac=NebFac, kcount_usr=kcount_usr, ExtremesFigSize=ExtremesFigSize, ResizeExtremesFac_max=ResizeExtremesFac_max, ResizeExtremesFac_min = ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)    



# Documentation: 
# "B" = array, initial guess for flux of "secondary"
# vrads1, vrads2 = RVs of primary, secondary
# vradsOut, vradsBin = RVs of tertiary * reflex motion of binary
# waves: wavelength grid on which disentanglement should take place
# Resid --> returns array of residual spectra between obs and dis1+dis2
# Reduce --> Returns reduced chi2
def disentangle3Comp(B, vrads1, vrads2, vradsOut, vradsBin, waves, ObsSpecs, weights, StrictNeg, PosLimCond, Poslimall,  nusdata, nusdataOut, Orbital_Params, K1s, K2s, KOuts, MJDs, phis, phisOut, specnames,  Rangestr, StarName, ScalingNeb, NebSpec, Resid=False, Reduce=False, ShowItr=False, Once=False,  InterKind='linear', itrnumlim=100, PLOTCONV=False, PLOTITR=False,   PLOTEXTREMES=False, PLOTFITS=False, kcount_extremeplot=0, linewidExt=3,  N_Iteration_Plot=50, NebLines = False, NebFac=1, kcount_usr=0, ExtremesFigSize=(7,8), ResizeExtremesFac_max=1., ResizeExtremesFac_min=1, TitleFileName=False, DisLine=None):    
    global kcount, k1, k2, kOut, DoFs, kcount, K1now, K2now, K3now
    C = waves*0.
    StrictNegA, StrictNegB, StrictNegC, StrictNegD = StrictNeg
# If convergence plot: allow spectra to be positive for sensible convergence plot:
    if PLOTCONV:
        StrictNegA, StrictNegB, StrictNegC, StrictNegD = False, False, False, False    
    PosLimCondA, PosLimCondB, PosLimCondC, PosLimCondD, PosLimCondNeb = PosLimCond
    if (not Once):
        try:
            K1now = K1s[k1]
            K2now = K2s[k2]   
            K3now = KOuts[kOut]   
            print("Disentangeling..... K1, K2, K3=", K1now, K2now, K3now)
        except:
            pass
    else:
        k1, k2, kOut = 0, 0, 0
    Facshift1 = np.sqrt( (1 + (vrads1+vradsBin)/clight) / (1 - (vrads1+vradsBin)/clight))
    Facshift2 = np.sqrt( (1 + (vrads2+vradsBin)/clight) / (1 - (vrads2+vradsBin)/clight))
    FacshiftOut = np.sqrt( (1 + vradsOut/clight) / (1 - vradsOut/clight))
    Ss1 = [interp1d(ObsSpecs[i][:, 0] /Facshift1[i], ObsSpecs[i][:, 1]-1.,
                      bounds_error=False, fill_value=0., kind=InterKind)(waves) for i in np.arange(len(ObsSpecs))] 
    Ss2 = [interp1d(ObsSpecs[i][:, 0] /Facshift2[i], ObsSpecs[i][:, 1]-1.,
                      bounds_error=False, fill_value=0., kind=InterKind)(waves) for i in np.arange(len(ObsSpecs))]      
    SsOut = [interp1d(ObsSpecs[i][:, 0] /FacshiftOut[i], ObsSpecs[i][:, 1]-1.,
                      bounds_error=False, fill_value=0., kind=InterKind)(waves) for i in np.arange(len(ObsSpecs))]      
#Frame of Refernce star 1:    
    WavesBA = np.array([waves * Facshift2[i] /Facshift1[i] for i in np.arange(len(vrads1))])  
    WavesCA = np.array([waves * FacshiftOut[i] /Facshift1[i] for i in np.arange(len(vrads1))])  
#Frame of Refernce star 2:        
    WavesAB = np.array([waves * Facshift1[i] /Facshift2[i] for i in np.arange(len(vrads1))])  
    WavesCB = np.array([waves * FacshiftOut[i] /Facshift2[i] for i in np.arange(len(vrads1))])  
#Frame of Refernce star 3:        
    WavesAC = np.array([waves * Facshift1[i] /FacshiftOut[i] for i in np.arange(len(vrads1))])  
    WavesBC = np.array([waves * Facshift2[i] /FacshiftOut[i] for i in np.arange(len(vrads1))])  
    if NebLines:
        FacshiftNeb = np.ones(len(vrads1))
        SsNeb = [interp1d(ObsSpecs[i][:, 0] / FacshiftNeb[i], ObsSpecs[i][:, 1]-1.,
                      bounds_error=False, fill_value=0., kind=InterKind)(waves) for i in np.arange(len(ObsSpecs))]             
        WavesNebA = np.array([waves *FacshiftNeb[i]/Facshift1[i] for i in np.arange(len(vrads1))])  
        WavesNebB = np.array([waves * FacshiftNeb[i] /Facshift2[i] for i in np.arange(len(vrads1))])         
        WavesNebC = np.array([waves * FacshiftNeb[i] /FacshiftOut[i] for i in np.arange(len(vrads1))])         
        WavesANeb = np.array([waves *Facshift1[i] /FacshiftNeb[i] for i in np.arange(len(vrads1))])  
        WavesBNeb = np.array([waves * Facshift2[i] /FacshiftNeb[i] for i in np.arange(len(vrads1))])   
        WavesCNeb = np.array([waves * FacshiftOut[i] /FacshiftNeb[i] for i in np.arange(len(vrads1))])   
    itr = 0
    while itr<itrnumlim:
        itr+=1
        BAshifts = [interp1d(WavesBA[i], B,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]    
        CAshifts = [interp1d(WavesCA[i], C,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]    
        if NebLines:
            NebAshifts = [interp1d(WavesNebA[i], NebSpec,bounds_error=False, fill_value=0., kind=InterKind) for i in np.arange(len(ObsSpecs))]    
            SpecMean = np.sum(np.array([[i]*(Ss1[i] - BAshifts[i](waves) - CAshifts[i](waves) - NebFac*ScalingNeb[i]*NebAshifts[i](waves)) for i in np.arange(len(Ss1))]), axis=0)      
        else:
            SpecMean = np.sum(np.array([weights[i]*(Ss1[i] - BAshifts[i](waves) - CAshifts[i](waves) ) for i in np.arange(len(Ss1))]), axis=0)          
        Anew = interp1d(waves, SpecMean, bounds_error=False, fill_value=0., kind=InterKind)(waves)  
        if StrictNegA:
            Anew = Limit(waves, Anew, Poslimall, PosLimCondA)        
        if 'A' in locals():
            Epsnew = np.amax((A - Anew)**2)        
        else:
            Epsnew = 0.
        A = np.copy(Anew)
        ABshifts = [interp1d(WavesAB[i], A,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]     
        CBshifts = [interp1d(WavesCB[i], C,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]     
        if NebLines:
            NebBshifts = [interp1d(WavesNebB[i], NebSpec,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]     
            SpecMean = np.sum(np.array([weights[i]*(Ss2[i] - ABshifts[i](waves) - CBshifts[i](waves) -  NebFac*ScalingNeb[i]*NebBshifts[i](waves)) for i in np.arange(len(Ss1))]), axis=0)
        else:
            SpecMean = np.sum(np.array([weights[i]*(Ss2[i] - ABshifts[i](waves) - CBshifts[i](waves)) for i in np.arange(len(Ss1))]), axis=0)                    
        Bnew = interp1d(waves, SpecMean, bounds_error=False, fill_value=0., kind=InterKind)(waves) 
        if StrictNegB:
            Bnew = Limit(waves, Bnew, Poslimall, PosLimCondB)            
        Epsnew = max(Epsnew, np.sum((B - Bnew)**2))
        B = np.copy(Bnew)
        ACshifts = [interp1d(WavesAC[i], A,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]     
        BCshifts = [interp1d(WavesBC[i], B,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]     
        if NebLines:
            NebCshifts = [interp1d(WavesNebC[i], NebSpec,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]     
            SpecMean = np.sum(np.array([weights[i]*(SsOut[i] - ACshifts[i](waves) - BCshifts[i](waves) -  NebFac*ScalingNeb[i]*NebCshifts[i](waves)) for i in np.arange(len(Ss1))]), axis=0)
        else:
            SpecMean = np.sum(np.array([weights[i]*(SsOut[i] - ACshifts[i](waves) - BCshifts[i](waves) ) for i in np.arange(len(Ss1))]), axis=0)                    
        Cnew = interp1d(waves, SpecMean, bounds_error=False, fill_value=0., kind=InterKind)(waves) 
        if StrictNegC:
            Cnew = Limit(waves, Cnew, Poslimall, PosLimCondC)            
        Epsnew = max(Epsnew, np.sum((C - Cnew)**2))
        C = np.copy(Cnew)       
        if NebLines:
            ANebshifts = [interp1d(WavesANeb[i], A,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]     
            BNebshifts = [interp1d(WavesBNeb[i], B,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]    
            CNebshifts = [interp1d(WavesCNeb[i], C,bounds_error=False, fill_value=0., kind=InterKind)
                        for i in np.arange(len(ObsSpecs))]                 
            SpecMean = np.sum(np.array([weights[i]*Limit(waves, SsNeb[i] - ANebshifts[i](waves) - BNebshifts[i](waves) - CNebshifts[i](waves), Poslimall, PosLimCondNeb) for i in np.arange(len(Ss1))]), axis=0)   
            NebSpecnew = interp1d(waves, SpecMean, bounds_error=False, fill_value=0., kind=InterKind)(waves)
            NebSpecnew[NebSpecnew<Poslimall[0]] = 0.
            Epsnew = max(Epsnew, np.sum(np.abs(NebSpec - NebSpecnew))) 
            NebSpec = NebSpecnew               
        if ShowItr:
            if (itr/itrnumlim*100)%10==0:
                print("Finished " + str(itr) + " out of " + str(itrnumlim) + " iterations (" +str(np.round(itr/itrnumlim*100.,3)) + "%)")
        if PLOTCONV:
            plt.scatter(itr, np.log10(Epsnew), color='blue')
        if PLOTITR:
            if itr%N_Iteration_Plot==0:
                plt.plot(waves, A, label=itr)
                plt.plot(waves, B, label=itr)
                plt.plot(waves, C, label=itr)
                if NebLines:
                    plt.plot(waves, NebSpec, label=itr)
    print("Finished after ", itr, " iterations")
    if PLOTCONV or PLOTITR:
        if PLOTITR:
            plt.ylabel('Normalised flux')
            plt.xlabel('Wavelength')
            plt.legend()
        if PLOTCONV:
            plt.ylabel('log(Eps)')
            plt.xlabel('iteration number')
        plt.show()  
    DisSpecVector = np.array([A, B , C, NebSpec])
    return DisSpecVector+1.,  CalcDiffs3D(DisSpecVector, vrads1,  vrads2,  vradsOut, vradsBin, waves, ObsSpecs, nusdata, nusdataOut, Orbital_Params, K1s, K2s, KOuts, MJDs, phis, phisOut, specnames, Rangestr, StarName, ScalingNeb, Resid=Resid, Reduce=Reduce, ShowItr=ShowItr,   PLOTEXTREMES=PLOTEXTREMES, PLOTFITS=PLOTFITS, kcount_extremeplot=kcount_extremeplot,  linewidExt=linewidExt, NebLines = NebLines, NebFac=NebFac, kcount_usr=kcount_usr, ExtremesFigSize=ExtremesFigSize, ResizeExtremesFac_max=ResizeExtremesFac_max, ResizeExtremesFac_min = ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)    

    
# Assuming spectrum for secondary (Bini) and vrads1, gamma, and K1, explore K2s array for best-fitting K2
# Ini = determines initial assumption for 0'th iteration
# ShowItr = determines whether 
def Grid_disentangling2D(waveRanges, nusdata, Bini, Orbital_Params, K1s, K2s,  ObsSpecs, weights, StrictNeg, PosLimCond, Poslimall,  MJDs, phis, specnames,  Rangestr, StarName, ScalingNeb, Ini=None, ShowItr=False,  InterKind='linear', itrnumlim = 100, NebOff=True, PLOTCONV=False, PLOTITR=False,   PLOTEXTREMES=False, PLOTFITS=False, kcount_extremeplot=0, linewidExt=3,ParbSize=3, N_Iteration_Plot=50,  NebLines = False, NebFac=1, kcount_usr=0, ExtremesFigSize=(7,8), ResizeExtremesFac_max=1., ResizeExtremesFac_min=1., TitleFileName=False, DisLine=None):
    global kcount, k1, k2, DoFs
    N = 0
    Diffs=np.zeros(len(K1s)*len(K2s)).reshape(len(K1s), len(K2s))  
    DoFs=0 
    for waves in waveRanges:       
        kcount = 0
        for k1, K1 in enumerate(K1s):
            for k2, K2 in enumerate(K2s):
                Orbital_Params_Updated = Orbital_Params.copy()
                Orbital_Params_Updated['K1'] = K1
                Orbital_Params_Updated['K2'] = K2
                vrads1, vrads2 = v1andv2(nusdata, Orbital_Params_Updated)  
                NebSpec = waves * 0.
                if Ini=='A':
                    print("Initial guess provided for component " + Ini)      
                    #print Bini(waves)
                    Diffs[k1,k2] += disentangle(Bini(waves), vrads2, vrads1,  waves, ObsSpecs, weights, StrictNeg, PosLimCond, Poslimall, nusdata, Orbital_Params_Updated, K1s, K2s, MJDs, phis, specnames, Rangestr, StarName,  ScalingNeb, NebSpec, InterKind=InterKind, itrnumlim = itrnumlim,  PLOTCONV=PLOTCONV, PLOTITR=PLOTITR,   PLOTEXTREMES=PLOTEXTREMES, PLOTFITS=PLOTFITS, kcount_extremeplot=kcount_extremeplot, linewidExt=linewidExt,  N_Iteration_Plot=N_Iteration_Plot, NebLines = NebLines, NebFac=1, kcount_usr=kcount_usr, ExtremesFigSize=ExtremesFigSize, ResizeExtremesFac_max=ResizeExtremesFac_max, ResizeExtremesFac_min = ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)[-1]      
                elif Ini=='B':
                    print("Initial guess provided for component " + Ini)        
                    Diffs[k1,k2] += disentangle(Bini(waves), vrads1, vrads2,  waves, ObsSpecs, weights, StrictNeg, PosLimCond, Poslimall, nusdata, Orbital_Params_Updated, K1s, K2s, MJDs, phis, specnames, Rangestr, StarName,  ScalingNeb, NebSpec, InterKind=InterKind, itrnumlim = itrnumlim,  PLOTCONV=PLOTCONV, PLOTITR=PLOTITR,   PLOTEXTREMES=PLOTEXTREMES, PLOTFITS=PLOTFITS, kcount_extremeplot=kcount_extremeplot, linewidExt=linewidExt, N_Iteration_Plot=N_Iteration_Plot, NebLines = NebLines, NebFac=1, kcount_usr=kcount_usr, ExtremesFigSize=ExtremesFigSize, ResizeExtremesFac_max=ResizeExtremesFac_max, ResizeExtremesFac_min = ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)[-1]              
                else:
                    print("No initial approximation given, assuming flat spectrum for secondary...")
                    Bini = interp1d(waves, np.ones(len(waves)),bounds_error=False, fill_value=1.)  
                    Diffs[k1,k2] += disentangle(Bini(waves), vrads2, vrads1,  waves, ObsSpecs, weights, StrictNeg, PosLimCond, Poslimall, nusdata, Orbital_Params_Updated, K1s, K2s, MJDs, phis, specnames, Rangestr, StarName,  ScalingNeb, NebSpec, InterKind=InterKind, itrnumlim = itrnumlim,  PLOTCONV=PLOTCONV, PLOTITR=PLOTITR,   PLOTEXTREMES=PLOTEXTREMES, PLOTFITS=PLOTFITS, kcount_extremeplot=kcount_extremeplot, linewidExt=linewidExt,  N_Iteration_Plot=N_Iteration_Plot, NebLines = NebLines, NebFac=1, kcount_usr=kcount_usr, ExtremesFigSize=ExtremesFigSize, ResizeExtremesFac_max=ResizeExtremesFac_max, ResizeExtremesFac_min=ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)[-1]         
    Diffs /= (DoFs)  
    try:
        StepSize1 = K1s[1] - K1s[0]
    except:
        StepSize1 = 0
    try:
        StepSize2 = K2s[1] - K2s[0]
    except:
        StepSize2 = 0        
    np.savetxt('disentangled/' + Rangestr + '_' + 'grid_dis_K1K2.txt', np.array(Diffs), header='#K1min, K2min, stepK1, stepK2, DoF = ' + str(K1s[0]) + ', ' + str(K2s[0]) + ', ' + str(StepSize1) + ', ' + str(StepSize2)   + ', ' + str(DoFs) ) 
    k1min, k2min = np.argwhere(Diffs == np.min(Diffs))[0]
    #print Diffs
    print("True velocities: ", k1min, k2min, K1s[k1min], K2s[k2min])
# Start with uncertainty on K1:    
    if StepSize2 > 0:
        chi2P, K2, K2err = Chi2con(Diffs[k1min,:], DoFs, K2s, Rangestr, comp='secondary', ParbSize=ParbSize)
    else:
        K2 = K2s[0]
        K2err = 0
    if StepSize1 > 0:
        chi2P, K1, K1err = Chi2con(Diffs[:,k2min], DoFs, K1s, Rangestr, comp='primary', ParbSize=ParbSize)
    else:
        K1 = K1s[0]     
        K1err = 0
    print("K1, K1 min error:", K1, K1err)         
    print("K2, K2 min error:", K2, K2err) 
    return K1, K2



    
# Assuming spectrum for secondary (Bini) and vrads1, gamma, and K1, explore K2s array for best-fitting K2
# Ini = determines initial assumption for 0'th iteration
# ShowItr = determines whether 
def Grid_disentangling3D(waveRanges, nusdata, nusdataOut, Bini, Orbital_Params, K1s, K2s,  KOuts, ObsSpecs, weights, StrictNeg, PosLimCond, Poslimall,  MJDs, phis, phisOut, specnames,  Rangestr, StarName, ScalingNeb, Ini=None, ShowItr=False,  InterKind='linear', itrnumlim = 100, NebOff=True, PLOTCONV=False, PLOTITR=False,   PLOTEXTREMES=False, PLOTFITS=False, kcount_extremeplot=0, linewidExt=3,  ParbSize=3, N_Iteration_Plot=50,  NebLines = False, NebFac=1, kcount_usr=0, ExtremesFigSize=(7,8), ResizeExtremesFac_max=1., ResizeExtremesFac_min=1., TitleFileName=False, DisLine=None):
    global kcount, k1, k2, kOut, DoFs
    N = 0
    Diffs=np.zeros(len(K1s)*len(K2s)*len(KOuts)).reshape(len(K1s), len(K2s), len(KOuts))  
    DoFs=0 
    for waves in waveRanges:       
        kcount = 0
        for k1, K1 in enumerate(K1s):
            for k2, K2 in enumerate(K2s):
                for kOut, KOut in enumerate(KOuts):
                    Orbital_Params_Updated = Orbital_Params.copy()
                    Orbital_Params_Updated['K1'] = K1
                    Orbital_Params_Updated['K2'] = K2
                    Orbital_Params_Updated['KOut'] = KOut
                    try:
                        vrads1, vrads2 = v1andv2(nusdata, Orbital_Params_Updated)  
                        vradsOut, vradsBin = v1andv2(nusdataOut, Orbital_Params_Updated, Config='Tertiary')  
                    except:
                        Diffs[k1,k2,kOut] = 1E5
                        kcount+=1
                        continue
                    NebSpec = waves * 0.
                    if Ini=='A':
                        print("Initial guess provided for component " + Ini)      
                        #print Bini(waves)
                        Diffs[k1,k2, kOut] += disentangle3Comp(Bini(waves), vrads2, vrads1, vradsOut, vradsBin,  waves, ObsSpecs, weights, StrictNeg, PosLimCond, Poslimall, nusdata, nusdataOut, Orbital_Params_Updated, K1s, K2s, KOuts, MJDs, phis, phisOut, specnames, Rangestr, StarName,  ScalingNeb, NebSpec, InterKind=InterKind, itrnumlim = itrnumlim,  PLOTCONV=PLOTCONV, PLOTITR=PLOTITR,   PLOTEXTREMES=PLOTEXTREMES, PLOTFITS=PLOTFITS, kcount_extremeplot=kcount_extremeplot, linewidExt=linewidExt,  N_Iteration_Plot=N_Iteration_Plot, NebLines = NebLines, NebFac=1, kcount_usr=kcount_usr, ExtremesFigSize=ExtremesFigSize, ResizeExtremesFac_max=ResizeExtremesFac_max, ResizeExtremesFac_min = ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)[-1]      
                    elif Ini=='B':
                        print("Initial guess provided for component " + Ini)        
                        Diffs[k1,k2, kOut] += disentangle3Comp(Bini(waves), vrads1, vrads2, vradsOut, vradsBin,   waves, ObsSpecs, weights, StrictNeg, PosLimCond, Poslimall, nusdata, nusdataOut, Orbital_Params_Updated, K1s, K2s, KOuts, MJDs, phis, phisOut, specnames, Rangestr, StarName,  ScalingNeb, NebSpec, InterKind=InterKind, itrnumlim = itrnumlim,  PLOTCONV=PLOTCONV, PLOTITR=PLOTITR,   PLOTEXTREMES=PLOTEXTREMES, PLOTFITS=PLOTFITS, kcount_extremeplot=kcount_extremeplot, linewidExt=linewidExt,  N_Iteration_Plot=N_Iteration_Plot, NebLines = NebLines, NebFac=1, kcount_usr=kcount_usr, ExtremesFigSize=ExtremesFigSize, ResizeExtremesFac_max=ResizeExtremesFac_max, ResizeExtremesFac_min = ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)[-1]                        
                    else:
                        print("No initial approximation given, assuming flat spectrum for secondary...")
                        Bini = interp1d(waves, np.ones(len(waves)),bounds_error=False, fill_value=1.)  
                        Diffs[k1,k2] += disentangle3Comp(Bini(waves), vrads2, vrads1, vradsOut, vradsBin,  waves, ObsSpecs, weights, StrictNeg, PosLimCond, Poslimall, nusdata, nusdataOut, Orbital_Params_Updated, K1s, K2s, KOuts, MJDs, phis, phisOut, specnames, Rangestr, StarName,  ScalingNeb, NebSpec, InterKind=InterKind, itrnumlim = itrnumlim,  PLOTCONV=PLOTCONV, PLOTITR=PLOTITR,   PLOTEXTREMES=PLOTEXTREMES, PLOTFITS=PLOTFITS, kcount_extremeplot=kcount_extremeplot, linewidExt=linewidExt, N_Iteration_Plot=N_Iteration_Plot, NebLines = NebLines, NebFac=1, kcount_usr=kcount_usr, ExtremesFigSize=ExtremesFigSize, ResizeExtremesFac_max=ResizeExtremesFac_max, ResizeExtremesFac_min=ResizeExtremesFac_min, TitleFileName=TitleFileName, DisLine=DisLine)[-1]         
    Diffs /= (DoFs)  
    try:
        StepSize1 = K1s[1] - K1s[0]
    except:
        StepSize1 = 0
    try:
        StepSize2 = K2s[1] - K2s[0]
    except:
        StepSize2 = 0        
    try:
        StepSize3 = KOuts[1] - KOuts[0]
    except:
        StepSize3 = 0        
    PickleHeader = {'KsMin': [K1s[0], K2s[0], KOuts[0]],
                    'Steps': [StepSize1, StepSize2, StepSize3],
                    'DoF' : DoFs }
    Chi2FileName = 'disentangled/' + Rangestr + '_' + 'grid_dis_K1K2K3.pkl'
    pickle.dump(PickleHeader, open(Chi2FileName, 'wb'))
    pickle.dump(Diffs, open(Chi2FileName, 'wb'))
    #np.savetxt('disentangled/' + Rangestr + '_' + 'grid_dis_K1K2K3.txt', (np.array(Diffs))., header='#K1min, K2min, K3min, stepK1, stepK2, stepK3, DoF = ' + str(K1s[0]) + ', ' + str(K2s[0]) + ', ' + ', ' + str(KOuts[0]) +  ', ' + str(StepSize1) + ', ' + str(StepSize2)   + ', ' + str(StepSize3) + ', ' + str(DoFs) )     
    k1min, k2min, kOutmin = np.argwhere(Diffs == np.min(Diffs))[0]
    #print Diffs
    print("True velocities: ", k1min, k2min, kOutmin, K1s[k1min], K2s[k2min], KOuts[kOutmin] )
# Start with uncertainty on K1:    
    if StepSize2 > 0:
        chi2P, K2, K2err = Chi2con(Diffs[k1min,:, kOutmin], DoFs,  K2s, Rangestr, comp='secondary', ParbSize=ParbSize)
    else:
        K2 = K2s[0]
        K2err = 0
    if StepSize1 > 0:
        chi2P, K1, K1err = Chi2con(Diffs[:,k2min, kOutmin], DoFs, K1s, Rangestr, comp='primary', ParbSize=ParbSize)
    else:
        K1 = K1s[0]     
        K1err = 0
    if StepSize3 > 0:
        chi2P, KOut, KOuterr = Chi2con(Diffs[k1min,k2min, :], DoFs, KOuts, Rangestr, comp='tertiary', ParbSize=ParbSize)
    else:
        KOut = KOuts[0]     
        KOuterr = 0
    print("K1, K1 min error:", K1, K1err)         
    print("K2, K2 min error:", K2, K2err) 
    print("KOut, KOut min error:", KOut, KOuterr) 
    return K1, K2, KOut
    
    

def read_FEROS(infile):
    print("%s: Input file is a FEROS file." % infile)
    header = fits.getheader(infile)
    try:
        flux = fits.getdata(infile)
        crval = header['CRVAL1']
        crpix = header['CRPIX1']
        cdelt = header['CDELT1']

        wave = crval - (cdelt * crpix - cdelt) + np.arange(flux.shape[0]) * cdelt
    except:
        table = Table.read(infile)
        flux = table['FLUX'][0]   
        wave = table['LAMBDA'][0]   
    return wave, flux


def read_GIRAFFE(infile):
    print(("%s: Input file is a GIRAFFE file." % infile))
    header = fits.getheader(infile)
    data = fits.getdata(infile)
    wl0 = header['CRVAL1']  # Starting wl at CRPIX1
    delt = header['CDELT1']  # Stepwidth of wl
    pix = header['CRPIX1']  # Reference Pixel
    wave = wl0 - (delt * pix - delt) + np.arange(data.shape[0]) * delt
    table = Table.read(infile, hdu=1)
    flux = table['NORM_SKY_SUB_CR']
    return wave, flux


# Flexible read_file function, credit @ Julia Bodensteiner
def read_file(infile, SkipLines=0):
    ext = str(infile.split('.')[-1])
    if (ext == 'fits') or (ext ==  'fit'):
        wave, flux = read_fits(infile)

    elif (ext == 'gz'):
        wave, flux = read_tlusty(infile)

    elif (ext == 'dat' or ext == 'ascii' or ext == 'txt' or ext == 'nspec'):
        wave, flux = read_xytable(infile, SkipLines=SkipLines)

    elif (ext == 'tfits'):
        wave, flux = read_uvespop(infile)

    elif (ext == 'hfits'):
        wave, flux = read_hermes_normalized(infile)

    else:
        wave, flux = read_xytable(infile, SkipLines=SkipLines)
    return np.array([wave, np.nan_to_num(flux,1.)]).T

def read_fits(infile):
    print("%s: Input file is a fits file." % infile)

    header = fits.getheader(infile)

    if 'HIERARCH SPECTRUM EXTRACTION' in header:
        wave, flux = read_psfSpec(infile)

    elif 'INSTRUME' in header:
        ins = header['INSTRUME']
        if (ins == 'MUSE'):
            wave, flux = read_pampelMUSE(infile)

        elif (ins == 'HERMES'):
            wave, flux = read_HERMES(infile)

        elif (ins == 'FEROS'):
            wave, flux = read_FEROS(infile)
        elif (ins == 'XSHOOTER'):
            wave, flux = read_XSHOOTER(infile)
        elif (ins == 'UVES'):
            wave, flux = read_UVES(infile)
        elif (ins == 'UVES_STITCH'):
            wave, flux = read_UVES_STITCH(infile)                    
        elif (ins == 'GIRAFFE' and 'nLR' in infile):
            wave, flux = read_GIRAFFE(infile)     
        elif (ins == 'GIRAFFE'):
            #wave, flux = read_GIRAFFE(infile)               
            wave, flux = read_GIRAFFE2(infile)               
        elif (ins == 'ESPCOUDE'):
            wave, flux = read_NLA(infile)   
        elif (ins == 'COS'):
            wave, flux = read_COS(infile)               
        elif (ins == 'STIS'):
            wave, flux = read_STIS(infile)              
        else:
            print('File type unkown, trying HERMES')
            wave, flux = read_HERMES(infile)                        
    else:
        wave, flux = read_HERMES(infile)

    return wave, flux

def read_xytable(infile, SkipLines=0):
    print(("%s: Input file is an xytable file." % infile))    
    #spec = (pd.read_csv(infile, sep=" ", header=0)).values   
    spec = (pd.read_csv(infile, header=None, delim_whitespace=True, skiprows=SkipLines)).values    
    wave = spec[:,0]
    flux = spec[:,1]    
    return wave, flux

def read_GIRAFFE2(infile):
    print("%s: Input file is a GIRAFFE2 file." % infile)
    hdunum=1
    header = fits.getheader(infile)
    hdul = fits.open(infile)
    table = Table.read(infile, hdu=hdunum)
    wave = table['WAVE'][0]*10.
    flux = table['FLUX_REDUCED'][0]    
    return wave, flux


def read_UVES_STITCH(infile):
    print("%s: Input file is a stitched UVES file" % infile)
    table = Table.read(infile, hdu=1)
    wave = table['WAVE']
    flux = table['FLUX']
    return wave, flux
        

def read_UVES(infile):
    print("%s: Input file is a UVES or GIRAFFE file." % infile)
    header = fits.getheader(infile)
    if 'SPEC_COM' in header:
        print('a')
        table = Table.read(infile, hdu=1)
        wave = table['wave']
        flux = table['FLUX_REDUCED']
    elif 'CRVAL1' in header:
        print('b')
        flux = fits.getdata(infile)
        crval = header['CRVAL1']
        cdelt = header['CDELT1']
        naxis1 = header['NAXIS1']
        wave = crval + np.arange(0, naxis1) * cdelt
    else:
        print('c')
        table = Table.read(infile, hdu=0)
        try:
            wave = table['WAVE'][0]
        except:
            wave = table['LAMBDA'][0] 
        try:
            flux = table['FLUX_REDUCED'][0]
        except:
            flux = table['FLUX'][0]    
    if wave[0] < 1000:  # wavelength is given in nm
        wave = wave * 10 # convert wavelength to A
    return wave, flux        
        

def read_pampelMUSE(infile):
    print(("%s: Input file is a pampelMUSE file." % infile))
    header = fits.getheader(infile)
    flux = fits.getdata(infile, ext=0)
    #err = fits.getdata(infile, ext=1)
    wl0 = header['CRVAL1']  # Starting wl at CRPIX1
    delt = header['CDELT1']  # Stepwidth of wl
    pix = header['CRPIX1']  # Reference Pixel
    wave = wl0 - (delt * pix - delt) + np.arange(flux.shape[0]) * delt

    # because pampelmuse gives the wl in m --> convert to A
    #wave = wave 
    #flux_out = flux, err

    return wave, flux


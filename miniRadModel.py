#Homebrew radiation code
# adapted heavily from Pierrehumbert's mini climate model

#This version includes  computation of
#the heating and flux due to stellar absorption.  The
#thermal infrared calculation has not changed.

#Notes:
#      --The temperature scaling for the continua assumes all continua
#        are defined at a reference temperature of 300K; this introduces
#        a slight error for the H2O continuum, which is defined relative
#        to 296K, but it's small compared to other issues with the
#        continuum. It should be fixed, eventually, though.
#


#Note regarding the stellar absorption: For simplicity
#the zenith angle is left at the same mean angle used
#in the thermal calculation.  It ought to be separated
#out as a separate variable, so that one can average
#over time if desired. In other words, instead of using
#the same cosThetaBar for thermal and stellar, this should
#be made into an argument of the stellar transmission function.
#
#Note also that this code does
#not include the effects of Rayleigh scattering.

#This is an elaborated version of miniClimt.py,
#which incorporates temperature weighting and allows
#for an inhomogeneous path, and for distinguishing
#self from foreign broadening (needed for continuum).
#It has a very approximate method for dealing with the
#distinction beween the self and foreign line broadening,
#incorporated through a single band-independent absorption
#ratio, called SelfRat.  This could easily be improved upon
#by reading in exponential sum tables for self and air broadening,
#and interpolating between the two (or using a band-dependent
#selfrat; if doing that, one might as well put in a band-dependent
#temperature scaling at the same time. For now, we've tried to keep
#things simple.
#
#Make sure to set selfrat consistently with the band table you
#are using. If you read in a band table that already was computed
#using self-broadening, then SelfRat should be set to 1 .
#
#This version doesn't include the semi-grey and Malkmus
#transmissions, which were in the simpler version
#only for pedagogical purposes.

#**ToDo:
#          *Speed up optically thin case
#          *Improve accuracy of optically thick case
#          *Put in temperature weighting (only here,not in
#              the simple version
#          *Write an alternate version based on
#           ODE integration (as prelude to scattering)
#          *Separate the functions out from the work loop,
#           so the radiation code can be used for multiple purposes
#          *The continuum as written assumes a self-continuum. This
#           is correct for water vapor, but not CO2. Need to find
#           a way to generalize that. Perhaps define separate self
#           and foreign continuum functions

#



import math
import string
import sys

#Things written for ClimateBook
sys.path.append('/home/jake/Documents/runaway/src/ppc/utilities/')
import phys
## needed just to read in the exponential sum tables, I think
from ClimateUtilities import *

#Mean slant path
cosThetaBar = .5

#Continuum functions (set continuum to desired function)
def NoContinuum(wave,T):
    return 0.
#
#CO2 foreign continuum
def CO2ForeignContinuum(wave,T):
    if (wave >= 25.) and (wave <= 550.):
        kappa = math.exp(
            -8.853 + 0.028534 *wave -0.00043194 *wave**2 + \
            1.4349e-6* wave**3  -1.5539e-9* wave**4)
    elif (wave >= 1050.) and (wave <= 1900.):
        kappa = math.exp(
        -537.09 +1.0886*wave -0.0007566 *wave**2 +1.8863e-7*wave**3\
        -8.2635e-12 *wave**4)
    else:
        kappa = 0.
    return (300./T)**1.7*kappa

#CO2 self continuum
def CO2SelfContinuum(wave,T):
    if (wave >= 25.) and (wave <= 550.):
        kappa = math.exp(
            -8.853 + 0.028534 *wave -0.00043194 *wave**2 + \
            1.4349e-6* wave**3  -1.5539e-9* wave**4)
    elif (wave >= 1050.) and (wave <= 1900.):
        kappa = math.exp(
        -537.09 +1.0886*wave -0.0007566 *wave**2 +1.8863e-7*wave**3\
        -8.2635e-12 *wave**4)
    else:
        kappa = 0.
    return 1.3*(300./T)**1.7*kappa

#Water vapor self-continuum
def H2OSelfContinuum(wave,T):
    if (wave>= 500.) and (wave<= 1400.):
        kappaC = math.exp(12.167-0.050898*wave +8.3207e-05*wave*wave
            -7.0748e-08*wave**3 + 2.3261e-11*wave**4)
    elif (wave>=2100.) and (wave <= 3000.):
        x = wave - 2500.
        kappaC = math.exp(-6.0055 + -0.0021363*x + 6.4723e-07*x**2 +
            -1.493e-08*x**3 + 2.5621e-11*x**4 + 7.328e-14*x**5)
    else:
        kappaC = 0.
    #Temperature dependence:
    return kappaC*(296./T)**4.25
H2OForeignContinuum = NoContinuum

ForeignContinuum = NoContinuum #Default
SelfContinuum = NoContinuum

#Planck functions in wavenumber space (cm**-1)
def Planck(wavenum,T):
    return 100.*math.pi*phys.c*phys.B(100.*wavenum*phys.c,T)

#Temperature derivative of Planck function
#Not currently used, but might be useful for implementing the
#optically thick limit to improve accuracy
#def dPlanck(wavenum,T):
#    nu = 100.*wavenum*phys.c
#    u = min(phys.h*nu/(phys.k*T),500.) #To prevent overflow
#    dB = (2.*phys.h*nu**3/phys.c**2)*(phys.h*nu/phys.k)*(1/T**2)*math.exp(u)/(math.exp(u)-1.)**2
#    return math.pi*100.*phys.c*dB

#Routine to read in band data from table
#by wavenumber.
def loadExpSumTable(fileName):
    c = readTable(fileName)
    bandList = []
    headers = c.listVariables()
    binHeaders = [h for h in headers if 'dH' in h]
    logKappaHeaders = [h for h in headers if 'logKappa' in h]
    nHeaders = len(binHeaders)
    for i in range(nHeaders):
        h = logKappaHeaders[i]
        h1 = binHeaders[i]
        waves = [string.atof(h.split('.')[1]),string.atof(h.split('.')[2])]
        bandParams = Dummy()
        bandParams.nu1 = string.atof(h.split('.')[1])
        bandParams.nu2 = string.atof(h.split('.')[2])
        bandParams.bandType = 'ExponentialSum'
        bandParams.kappa = Numeric.exp(c[h])
        bandParams.dH = c[h1]
        #
        bandList.append(bandParams)
    return bandList



#Band-averaged transmission function
#To keep track of self-collision vs. foreign collisions,
#we need two different pressure-weighted paths. "path"
#is the foreign-collision path, whereas "pathq" is the
#self-collision path. If we were getting really fancy, we
#might need more paths, since the temperature dependence
#of the continuum is different from that of the line absorption
def TransEsums(path,pathq,bandParams):
    '''compute band-averaged transmission function for a given path, pathq'''
    #Speed this up by dropping small terms from the sum.
    #It's easy to find which are smalll, since kappa is
    #sorted in increasing order.
    #
    #Handle the continuum (Reference Temperature curently hard-wired to 300K)
    #Temperature scaling is currently handled in the path calculation
    #Fix this, to generalize (use separate self and foreign contin fn)
    nu = (bandParams.nu1 + bandParams.nu2)/2.
    continOpacity = \
                  ForeignContinuum(nu,300.)*path + SelfContinuum(nu,300.)*pathq
    contin = math.exp(-continOpacity)
    #
    #The following statement cuts off the exponential sum
    #when the argument of the exponential gets too large. It
    #prevents underflow and also speeds things up in the optically
    #thick case
    pathTot = path + SelfRat*pathq #**This is where we handle the enhancement of
                           #self-collisions.  Make this band-dependent
                           #for greater accuracy.
    imax = Numeric.searchsorted(bandParams.kappa*pathTot,15.)
    return contin*sum(bandParams.dH[:imax] * Numeric.exp(-pathTot*bandParams.kappa[:imax]))

def TransGrey(path, pathq, kappa):
  '''compute transmission function for a grey gas'''
  path = path + pathq
  return math.exp(-path*kappa)

def TransSemiGrey(path, pathq, bandParams):
  '''compute transmission function for semi-grey gas (one kappa per band)'''
  path = path + pathq
  return math.exp(-path*bandParams.kappa)

def transmission(p,T, q, bandParams):
  '''compute the band-averaged transmission function between space and a given height'''
  n = len(p)
  Delta = bandParams.nu2-bandParams.nu1
  wave = (bandParams.nu2+bandParams.nu1)/2.
  #
  #Mass ratio for computing molar concentration
  #  (used to compute proportion of self-collisions)
  massrat = GHG.MolecularWeight/BackgroundGas.MolecularWeight
  selfRat = 1.
  trans_up = Numeric.zeros(n,Numeric.Float)
  trans_down = Numeric.zeros(n,Numeric.Float)

  #Upward
  for i in range(n):
      path = 0.
      pathq = 0.
      for j in range(i,n-1):

          moleCon = q[j]/(q[j] + (1-q[j])*massrat)

          wTS, wTF = path_scale(T[j], wave)
          #**Move this to midpoint for accuracy
          path = path + wTF*((1-moleCon)*p[j]/pref)*q[j]*abs(p[j+1]-p[j])/(g*cosThetaBar)
          #Compute the self-broadened pressure path
          pathq = pathq + wTS*moleCon*(p[j]/pref)*q[j]*abs(p[j+1]-p[j])/(g*cosThetaBar)
          if j == i:
            Trans1 = Trans(path,pathq,bandParams)
          #This will be inaccurate for j=i when an individual
          #layer becomes optically thick
          if Trans1 < 1.e-4:
              break
      trans_up[i] = Trans1 # trans at j == i (transmissivity from this layer to TOA)


  #Downward flux
  for i in range(n):
      path = 0.
      pathq= 0.
      for j in range(i,0,-1):
          #Note: we can use dB directly instead of
          #          dB/dT * dT/dp * dp
          moleCon = q[j]/(q[j] + (1-q[j])*massrat)
          wTS, wTF = path_scale(T[j], wave)
          #---------------------------------
          #**Fix this; move to midpoint
          path = path + wTF*((1-moleCon)*p[j]/pref)*q[j]*abs(p[j-1]-p[j])/(g*cosThetaBar)
          #Compute the self-broadened pressure path
          pathq = pathq + wTS*moleCon*(p[j]/pref)*q[j]*abs(p[j-1]-p[j])/(g*cosThetaBar)
          if j == i:
             Trans1 = Trans(path,pathq,bandParams)
          #Trans1 = Trans(path,pathq,bandParams)
          if Trans1 < 1.e-4:
              break
      trans_down[i] = Trans1 # trans at j = i (to surface?)

  return trans_up, trans_down

def path_scale(Tlev, wave):
  '''returns the self-broadened (wTF) and foreign broadened paths (wTs)'''
  #Temperature scaling------------
  #As currently written we
  #have one Tstar for all bands.  Further, if there
  #is a continuum, the scaling is taken from the continuum,
  #whereas it would be better to introduce a separate
  #path for the continuum
  if ForeignContinuum(wave,300.) > 1.e-6:
      wTF = ForeignContinuum(wave,Tlev)/ForeignContinuum(wave,300.)
  else:
      wTF = math.exp(-Tstar*(1./Tlev - 1./Tref))
  if SelfContinuum(wave,300.) > 1.e-6:
      wTS = SelfContinuum(wave,Tlev)/SelfContinuum(wave,300.)
  else:
      wTS = math.exp(-Tstar*(1./Tlev - 1./Tref))
  return wTS, wTF

#Note: q is now an array
#In this version, we compute the path in the transmission
#loop, which is more efficient than doing the whole integral
#from scratch in the transmission function routine.
#Note:All the path computations need to be cleaned up


def LWFluxBand(p,T,Tg,q,bandParams):
    '''returns LW flux_up, LW flux_down, and surface contribution to the OLR for a band'''
    n = len(p)
    Delta = bandParams.nu2-bandParams.nu1
    wave = (bandParams.nu2+bandParams.nu1)/2.
    #
    #Mass ratio for computing molar concentration
    #  (used to compute proportion of self-collisions)
    massrat = GHG.MolecularWeight/BackgroundGas.MolecularWeight
    #Temporary arrays for flux computation
    flux_up = Numeric.zeros(n,Numeric.Float)
    flux_down = Numeric.zeros(n,Numeric.Float)
    #-------------------------------------
    #Upward flux
    for i in range(n):
        Sum = 0.
        path = 0.
        pathq = 0.
        for j in range(i,n-1):
            #Note: we can use dB directly instead of
            #          dB/dT * dT/dp * dp
            dB = Planck(wave,T[j+1])-Planck(wave,T[j])
            moleCon = q[j]/(q[j] + (1-q[j])*massrat)

            wTS, wTF = path_scale(T[j], wave)
            #**Move this to midpoint for accuracy
            path = path + wTF*((1-moleCon)*p[j]/pref)*q[j]*abs(p[j+1]-p[j])/(g*cosThetaBar)
            #Compute the self-broadened pressure path
            pathq = pathq + wTS*moleCon*(p[j]/pref)*q[j]*abs(p[j+1]-p[j])/(g*cosThetaBar)
            Trans1 = Trans(path,pathq,bandParams)
            #This will be inaccurate for j=i when an individual
            #layer becomes optically thick
            Sum += Trans1*dB*Delta
            if Trans1 < 1.e-4:
                break

        #Add in boundary term
        BddTerm = (Planck(wave,Tg)-Planck(wave,T[-1]))*Trans(path,pathq,bandParams)
        flux_up[i] = Sum + (Planck(wave,T[i]) + BddTerm)*Delta
        # record Surface contribution to OLR at TOA
        if i == 0:
          SurfOLR = Planck(wave, Tg)*Trans(path, pathq, bandParams)*Delta

    #Downward flux
    for i in range(n):
        Sum = 0.
        path = 0.
        pathq= 0.
        for j in range(i,0,-1):
            #Note: we can use dB directly instead of
            #          dB/dT * dT/dp * dp
            dB = Planck(wave,T[j-1])- Planck(wave,T[j])
            moleCon = q[j]/(q[j] + (1-q[j])*massrat)
            wTS, wTF = path_scale(T[j], wave)
            #---------------------------------
            #**Fix this; move to midpoint
            path = path + wTF*((1-moleCon)*p[j]/pref)*q[j]*abs(p[j-1]-p[j])/(g*cosThetaBar)
            #Compute the self-broadened pressure path
            pathq = pathq + wTS*moleCon*(p[j]/pref)*q[j]*abs(p[j-1]-p[j])/(g*cosThetaBar)
            Trans1 = Trans(path,pathq,bandParams)
            Sum += Trans1*dB*Delta
            if Trans1 < 1.e-4:
                break
        #Add in boundary term
        BddTerm = (-Planck(wave,T[0]))*Trans(path,pathq,bandParams)
        flux_down[i] = Sum + (Planck(wave,T[i]) + BddTerm)*Delta
    return flux_up, flux_down, SurfOLR

def HeatBand(flux_up, flux_down, p):
    '''compute bandwise heating as gradient of fluxes using second-order
    centered differences,
    NOT verified'''
    n = len(p)
    heat = Numeric.zeros(n,Numeric.Float)
    flux = flux_up - flux_down
    #
    #Compute heat by taking the gradient of the flux.
    #The second order centered difference below allows
    #the computation to be done accurately even if
    #the pressure level spacing is non-uniform.
    #
    #Returns heating rate as W/kg.
    #Divide by specific heat to get K/s
    #
    delPlus = p[2:] - p[1:-1] #Should wind up as a [1:-1] array
    delMinus = p[1:-1]-p[:-2] #likewise
    A = (delMinus/delPlus)/(delPlus+delMinus)
    B = 1./delMinus - 1./delPlus
    C = (delPlus/delMinus)/(delPlus+delMinus)
    heat[1:-1] = A*flux[2:] + B*flux[1:-1] - C*flux[:-2]
    heat = g*heat #Convert to W/kg
    heat[0] = heat[1]
    heat[-1] = heat[-2]
    return heat

def SWFluxBand(p,T,Tg,q,bandParams):
    n = len(p)
    Delta = bandParams.nu2-bandParams.nu1
    wave = (bandParams.nu2+bandParams.nu1)/2.
    #
    #Mass ratio for computing molar concentration
    #  (used to compute proportion of self-collisions)
    massrat = GHG.MolecularWeight/BackgroundGas.MolecularWeight

    #Now compute the transmission function between each
    #level and the top.
    OrbitFact = (Lstellar/4.)/(phys.sigma*Tstellar**4)
    swFlux = OrbitFact*Delta*Planck(wave,Tstellar)
    flux = Numeric.zeros(n,Numeric.Float)
    flux[0] = -swFlux #Downward flux is negative
    #
    path = 0.
    pathq= 0.
    for j in range(1,n):
        moleCon = q[j]/(q[j] + (1-q[j])*massrat)
        #Temperature scaling------------
        #As currently written we
        #have one Tstar for all bands.  Further, if there
        #is a continuum, the scaling is taken from the continuum,
        #whereas it would be better to introduce a separate
        #path for the continuum
        if ForeignContinuum(wave,300.) > 1.e-6:
            wTF = ForeignContinuum(wave,T[j])/ForeignContinuum(wave,300.)
        else:
            wTF = math.exp(-Tstar*(1./T[j] - 1./Tref))
        if SelfContinuum(wave,300.) > 1.e-6:
            wTS = SelfContinuum(wave,T[j])/SelfContinuum(wave,300.)
        else:
            wTS = math.exp(-Tstar*(1./T[j] - 1./Tref))
        #---------------------------------
        #**Fix this; move to midpoint
        path = path + wTF*((1-moleCon)*p[j]/pref)*q[j]*abs(p[j-1]-p[j])/(g*cosThetaBar)
        #Compute the self-broadened pressure path
        pathq = pathq + wTS*moleCon*(p[j]/pref)*q[j]*abs(p[j-1]-p[j])/(g*cosThetaBar)
        Trans1 = Trans(path,pathq,bandParams)
        flux[j] = -swFlux*Trans1
    return flux

def LWtotal(p, T, Tg, q):
  '''return net flux and heating across all LW bands'''
  n = len(p)
  flux = Numeric.zeros(n,Numeric.Float)
  heat = Numeric.zeros(n,Numeric.Float)
  for band in bandData:
      if band.nu1 < LWCutoff:
          flux_up, flux_down, _ = LWFluxBand(p, T, Tg, q, band)
          heat_band = HeatBand(flux_up, flux_down, p)
          flux += (flux_up-flux_down)
          heat += heat_band
  return flux, heat

def LWgrey(p, T, Tg, q):
  '''LW for grey radiation scheme; transmission function must be set to TransGrey'''
  massrat = GHG.MolecularWeight/BackgroundGas.MolecularWeight
  n = len(p)
  #Mass ratio for computing molar concentration
  #  (used to compute proportion of self-collisions)
  #Temporary arrays for flux computation
  Ip = Numeric.zeros(n,Numeric.Float)
  Im = Numeric.zeros(n,Numeric.Float)
  #-------------------------------------
  #Upward flux
  for i in range(n):
      Sum = 0.
      path = 0.
      pathq = 0.
      for j in range(i,n-1):
          #Note: we can use dB directly instead of
          #          dB/dT * dT/dp * dp
          dB = phys.sigma*T[j+1]**4- phys.sigma*T[j]**4
          moleCon = q[j]/(q[j] + (1-q[j])*massrat)
          #Temperature scaling------------
          #As currently written we
          #have one Tstar for all bands.  Further, if there
          #is a continuum, the scaling is taken from the continuum,
          #whereas it would be better to introduce a separate
          #path for the continuum
          wTF = math.exp(-Tstar*(1./T[j] - 1./Tref))
          wTS = math.exp(-Tstar*(1./T[j] - 1./Tref))
          #---------------------------------
          #**Move this to midpoint for accuracy
          path = path + wTF*((1-moleCon)*p[j]/pref)*q[j]*abs(p[j+1]-p[j])/(g*cosThetaBar)
          #Compute the self-broadened pressure path
          pathq = pathq + wTS*moleCon*(p[j]/pref)*q[j]*abs(p[j+1]-p[j])/(g*cosThetaBar)
          Trans1 = Trans(path,pathq,kappa)
          #This will be inaccurate for j=i when an individual
          #layer becomes optically thick
          Sum += Trans1*dB
          if Trans1 < 1.e-4:
              break

      #Add in boundary term
      BddTerm = (phys.sigma*Tg**4 - phys.sigma*T[-1]**4)*Trans(path,pathq,kappa)
      #BddTerm = (Planck(wave,Tg)-Planck(wave,T[-1]))*Trans(path,pathq,bandParams)
      #SurfOLR = Planck(wave, Tg)*Trans(path, pathq, bandParams)
      Ip[i] = Sum + (phys.sigma*T[i]**4 + BddTerm)



  #Downward flux
  for i in range(n):
      Sum = 0.
      path = 0.
      pathq= 0.
      for j in range(i,0,-1):
          #Note: we can use dB directly instead of
          #          dB/dT * dT/dp * dp
          dB = phys.sigma*T[j-1]**4- phys.sigma*T[j]**4
          moleCon = q[j]/(q[j] + (1-q[j])*massrat)
          #Temperature scaling------------
          #As currently written we
          #have one Tstar for all bands.  Further, if there
          #is a continuum, the scaling is taken from the continuum,
          #whereas it would be better to introduce a separate
          #path for the continuum
          wTF = math.exp(-Tstar*(1./T[j] - 1./Tref))
          wTS = math.exp(-Tstar*(1./T[j] - 1./Tref))
          #---------------------------------
          #**Fix this; move to midpoint
          path = path + wTF*((1-moleCon)*p[j]/pref)*q[j]*abs(p[j-1]-p[j])/(g*cosThetaBar)
          #Compute the self-broadened pressure path
          pathq = pathq + wTS*moleCon*(p[j]/pref)*q[j]*abs(p[j-1]-p[j])/(g*cosThetaBar)
          Trans1 = Trans(path,pathq,kappa)
          Sum += Trans1*dB
          if Trans1 < 1.e-4:
              break
      #Add in boundary term
      BddTerm = (-phys.sigma*T[0]**4)*Trans(path,pathq,kappa)
      Im[i] = Sum + (phys.sigma*T[i]**4 + BddTerm)
  flux = Ip-Im
  heat = Numeric.zeros(n,Numeric.Float)
  #
  #Compute heat by taking the gradient of the flux.
  #The second order centered difference below allows
  #the computation to be done accurately even if
  #the pressure level spacing is non-uniform.
  #
  #Returns heating rate as W/kg.
  #Divide by specific heat to get K/s
  #
  delPlus = p[2:] - p[1:-1] #Should wind up as a [1:-1] array
  delMinus = p[1:-1]-p[:-2] #likewise
  A = (delMinus/delPlus)/(delPlus+delMinus)
  B = 1./delMinus - 1./delPlus
  C= (delPlus/delMinus)/(delPlus+delMinus)
  heat[1:-1] = A*flux[2:] + B*flux[1:-1] - C*flux[:-2]
  heat = g*heat #Convert to W/kg
  heat[0]=heat[1]
  heat[-1] = heat[-2]
  return flux, heat

## DEPRECATED; ADD GENERIC FUNCTION TO SUM OVER BANDS
#Sum of flux and heating over all bands
# def radCompLW(p,T,Tg,q):
#     n = len(p)
#     flux = Numeric.zeros(n,Numeric.Float)
#     heat = Numeric.zeros(n,Numeric.Float)
#     for band in bandData:
#         if band.nu1 < LWCutoff:
#             fluxBand,heatBand = LWFluxBand(p,T,Tg,q,band)
#             flux += fluxBand
#             heat += heatBand
#     return flux,heat
#
# def radCompStellar(p,T,Tg,q):
#     n = len(p)
#     flux = Numeric.zeros(n,Numeric.Float)
#     heat = Numeric.zeros(n,Numeric.Float)
#     for band in bandData:
#         #Shortwave (stellar) contribution
#         fluxBand,heatBand = SWFluxBand(p,T,Tg,q,band)
#         flux += fluxBand
#         heat += heatBand
#     return flux,heat

#----------------End of radiation model------------------
#
#Misc. utility functions
#
#OLR function for dry air adiabat
#Although our basic radiation calculation uses mass
#concentrations, this one takes molar concentrations
#of the greenhouse gas (in ppmv) as input, to make the
#units more familiar.  It also increases surface pressure
#to account for the additional mass of greenhouse gas.
#The surface pressure adjustment makes little difference out
#to 10% CO2, but as one goes to 20% and above it starts to become
#important
#
#Do we want to put in an isothermal stratosphere, to make
#this function more parallel to the ccmrad function?
def OLRDryAir(psAir,Tg,GHGmolarcon):
    ptop = 100.
    ps = psAir/(1.-1.e-6*GHGmolarcon)
    Mbar = 1.e-6*GHGmolarcon*GHG.MolecularWeight + (1.-1.e-6*GHGmolarcon)*BackgroundGas.MolecularWeight
    p = setpLin(ps,ptop,40)
    T = Tg*(p/ps)**(2./7.) #Ignores effect of GHG on R/cp
    q = 1.e-6*GHGmolarcon*GHG.MolecularWeight/Mbar
    q = q*Numeric.ones(40,Numeric.Float)
    flux,heat = radComp(p,T,Tg,q)
    #Add in contribution of part of spectrum outside the band table
    nu1 = bandData[0].nu1
    nu2 = bandData[-1].nu2
    m = romberg(Planck) #Evaluates definite integral
    flux1 = m([.001,nu1],Tg) + m([nu2,50000.],Tg)
    return flux[0] + flux1
#Functions for setting up pressure arrays
#
#Linear pressure array
def setpLin(ps,ptop,n):
    p = [ptop + i*(ps-ptop)/(n-1) for i in range(n)]
    return Numeric.array(p)

#Log pressure array (more resolution in stratosphere)
def setpLog(ps,ptop,n):
    rat = (ps/ptop)**(1./(n-1))
    p = [ptop*rat**i for i in range(n)]
    return Numeric.array(p)

#Linear array with extra resolution near ground
def setpExtraGroundRes(ps,ptop,n):
    n1 = 2*n/3
    n2 = n-n1
    p = [ptop + ((.9*ps-ptop)*i)/(n1-1) for i in range(n1)]
    p1 = [.9*ps + .1*ps*i/(n2) for i in range(1,n2+1)]
    p = p + p1 #Concatenation, not addition!
    return Numeric.array(p)

#-----------------------------------------------------------------------
#Remember to set gravity for your planet! It can be over-ridded
g = 9.8

#Read in band data (Moved to the calling program)
#


#Default data
BackgroundGas = phys.air #The transparent background gas
GHG =phys.CO2 #The greenhouse gas
#
pref = 1.e4 #Reference pressure at which band data is given
Tref = 260. #Reference temperature at which band data is given
            #Used for in-band scaling. Continuum has separate
            #temperature scaling.
            #**ToDo: Eliminate these defaults and force user to choose
Trans = TransEsums
SelfRat = 1. #Set mean ratio of self to foreign broadening for the
             #gas in use

#Set the continuum, if there is one
#SelfContinuum = CO2SelfContinuum
#ForeignContinuum = CO2ForeignContinuum

#To speed up the calculation, when doing both
#stellar absorption and longwave flux together,
#you can set a cutoff wavenumber for the longwave
#flux. That can be helpful, since the atmosphere
#is usually a lot cooler than the photosphere and
#therefore doesn't emit much at the shorter
#waves that are important for the incoming stellar
#absorption. Use with caution, though: if the
#atmosphere gets very hot, then you need the shorter
#wave emission
LWCutoff = 1.e6 #Wavenumber cutoff in 1/cm for longwave calc

#Set the photosphere temperature of the star
#(for shortwave radiation calculation
Tstellar = 3500.
#Set the stellar constant at the orbit
Lstellar = 4.*700.

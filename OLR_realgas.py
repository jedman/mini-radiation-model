
# TO DO
### generic band-sum function (takes fn that returns spectral, returns sum over wavenumber)
import numpy as np
import scipy.interpolate as sp
import matplotlib.pyplot as plt
import sys

home_dir = '/Users/jedman/Google Drive/Clouds/code/'
sys.path.append(home_dir +'runaway/src/')
#import miniClimtFancy as rad
import miniClimtFancy_jake as rad
import phys

#Path to the workbook datasets
datapath = '/Users/jedman/Documents/WorkbookDatasets/'
#Path to the exponential sum tables
EsumPath = datapath + 'Chapter4Data/ExpSumTables/'


#Initialize the radiation calculation, and over-ride defaults
#For doing the pure water vapor case, we read in the self-broadened
#table and leave rad.SelfRat = 1
#rad.bandData = rad.loadExpSumTable(EsumPath+'H2OTable.260K.100mb.self.data')
#
#For doing air/water mixtures, we read in the air-broadened table,
#and use a uniform rad.SelfRat. For high water vapor concentrations,
#this is less accurate than the previous method.
rad.bandData = rad.loadExpSumTable(EsumPath + 'H2OTable.260K.100mb.air.data')
rad.SelfRat = 5.

rad.ForeignContinuum = rad.NoContinuum #No w.v. foreign continuum
rad.BackgroundGas = phys.N2 #The transparent background gas
rad.GHG = phys.H2O #The greenhouse gas
rad.Tstar = 0. #Ignore temperature scaling except for continuum
rad.SelfContinuum = rad.H2OSelfContinuum

def OLR_spectral(Tg, psAir=1.e5, ptop = 0.01, n = 20, rh = 0.7):
  '''OLR as a function of wavenumber '''
  p, T, q = make_atmosphere(Tg, psAir = psAir, ptop = ptop, n = n, rh = rh)
  wave = []
  OLR = []
  surf = []
  for band in rad.bandData:
      flux_up, flux_down, surfOLR = rad.LWFluxBand(p,T,Tg,q,band)
      flux = flux_up - flux_down
      wave1 = (band.nu2+band.nu1)/2.
      Delta = (band.nu2-band.nu1)
      wave.append(wave1)
      #Convert flux at TOA to per wavenumber
      OLR.append(flux[0]/Delta)
      surf.append(surfOLR/Delta)
  return wave, OLR, surf

def OLR_spectral_atm(Tg, make_atm,  **kwargs):
  '''OLR as a function of wavenumber '''
  p, T, q = make_atm(Tg, **kwargs) # make
  wave = []
  OLR = []
  surf = []
  for band in rad.bandData:
      flux_up, flux_down, surfOLR = rad.LWFluxBand(p,T,Tg,q,band)
      flux = flux_up - flux_down
      wave1 = (band.nu2+band.nu1)/2.
      Delta = (band.nu2-band.nu1)
      wave.append(wave1)
      #Convert flux at TOA to per wavenumber
      OLR.append(flux[0]/Delta)
      surf.append(surfOLR/Delta)
  return wave, OLR, surf

def band_heating_atm(Tg, make_atm,  **kwargs):
  '''returns tuple of wavenumber, spectrally-resolved heating, p, T, q'''
  p, T, q = make_atm(Tg, **kwargs)
  fluxes = [rad.LWFluxBand(p, T, Tg, q, band) for band in rad.bandData]
  band_heating = np.array([rad.HeatBand(x[0], x[1], p) for x in fluxes])
  wave = np.array([(band.nu2 + band.nu1)/2 for band in rad.bandData])
  return wave, band_heating, p, T, q

def band_heating(Tg, n = 20, **kwargs):
  '''returns tuple of wavenumber, spectrally-resolved heating, p, T, q'''
  p, T, q = make_atmosphere(Tg, n = n, **kwargs)
  fluxes = [rad.LWFluxBand(p, T, Tg, q, band) for band in rad.bandData]
  band_heating = np.array([rad.HeatBand(x[0], x[1], p) for x in fluxes])
  wave = np.array([(band.nu2 + band.nu1)/2 for band in rad.bandData])
  return wave, band_heating, p, T, q


# utility
##for making a more realistic atmospheric water vapor profile
def make_atmosphere(Tg , psAir = 1.e5, ptop = 0.1, n = 20,  rh = 0.7, Ttrop = 200):
  ''' make an atmosphere from a moist adiabat
  TO DO: add stratosphere for temps greater than Ttrop'''
  psat = phys.satvps_function(phys.water)
  ps = psAir + psat(Tg) #Total surface pressure for this temperature
  p = rad.setpLog(ps,ptop,n)
  #p = rad.setpLin(ps,ptop,n)

  m = phys.MoistAdiabat()
  m.ptop = ptop
  p,T,molarCon, q = m(psAir,Tg,p)
  q[np.argwhere(T>Ttrop)] = rh*q[np.argwhere(T>Ttrop)]
  return p, T, q

def sum_over_bands(olr, wave):
  '''sums spectral olr over wavenumber'''
  return sum(olr)*(wave[1] - wave[0])


def make_atm_const(Tg, lapse_rate = 6. , psAir = 1.e5, ptop = 0.1, n = 20, rh = 0.7):
  '''make an atmosphere with a constant lapse rate (K/km)'''
  psat = phys.satvps_function(phys.water)
  ps = psAir + psat(Tg) #Total surface pressure for this temperature
  #ps = psAir # neglect psat(Tg) for now
  p = rad.setpLog(ps,ptop,n)
  z = np.zeros(n)
  rhosurf = 1.1455 # kg/m^3

  # make T(z), used for interpolating to pgrid from pierrehumbert
  ztop = 48.
  ggr = 9.81
  zlev = np.linspace(0,ztop, 1000)
  dz = zlev[1] - zlev[0] # (km)
  Tlev = np.zeros(zlev.shape)
  for i, z in enumerate(zlev):
      if i == 0:
          Tlev[i] = Tg
      else:
          Tlev[i] = Tlev[i-1] - dz*lapse_rate
  pqlev = np.zeros(zlev.shape)
  # make pq(zlev)
  pqlev = [psat(t) for t in Tlev]
  Ra = phys.air.R # J/kg K
  plev = [(ps)*(1 - lapse_rate*zi/Tg)**(ggr/Ra*1000/lapse_rate) for zi in zlev]
  qlev = rh*np.array(pqlev)/np.array(plev)
  #return plev, Tlev, qlev, zlev

  # interpolate to pressure levels
  # interpolate to pressure levels
  T = sp.interp1d(plev, Tlev, kind = 'linear', bounds_error = False, fill_value=(Tlev[-1], Tlev[0]))(p)
  q = sp.interp1d(plev, qlev,kind= 'linear', bounds_error = False, fill_value=(qlev[-1], qlev[0]))(p)
  return p, T, q

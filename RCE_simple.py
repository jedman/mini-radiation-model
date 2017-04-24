#import sys
import math
import numpy as np
import matplotlib.pyplot as plt
#sys.path.append('/Users/jedman/Google Drive/Clouds/code/runaway/src/')
import phys
import miniRadModel as rad

# semi-grey band data
# semiGreyBandParams = Dummy()
# semiGreyBandParams.nu1 = 650.
# semiGreyBandParams.nu2 = 700.
# semiGreyBandParams.kappa = 1.
# rad.bandData = [semiGreyBandParams]

#rad.bandData = rad.loadExpSumTable(EsumPath + 'H2OTable.260K.100mb.air.data')
rad.SelfRat = 1.

rad.ForeignContinuum = rad.NoContinuum #No w.v. foreign continuum
rad.BackgroundGas = phys.N2 #The transparent background gas
rad.GHG = phys.H2O #The greenhouse gas
rad.Tstar = 0. #Ignore temperature scaling
rad.SelfContinuum = rad.NoContinuum
rad.g = 9.81
rad.pref = 1.e4

n = 100 #n = 60
ps = 1.e5
ptop = 10. #Use 100 except for Venus, where we use 1.e4
#p = rad.setpExtraGroundRes(ps,ptop,n)
p = rad.setpLin(ps,ptop,n)
#p = rad.setpLog(ps,ptop,n) #Puts extra resolution in stratosphere.

class Dummy:
  pass

class RCE():
  '''do a semi-grey or grey rce'''
  def __init__(self, kappa = 1. , rh = 1., ps = 1.e5):
      self.kappa = kappa
      self.rh = rh
      self.dtime = 50*86400/1003 #  timestep and conversion to K/day
      self.massrat = rad.GHG.MolecularWeight/rad.BackgroundGas.MolecularWeight
      self.ps = 1.e5
      self.data = None
      self.psat = phys.satvps_function(phys.water)
      return

  def doRCE(self, Tglist, rce_type = 'semigrey'):
      '''calculate RCE for list of Tg'''
      if rce_type == 'grey':
          print "doing grey RCE calculation"
          #set grey kappa
          rad.kappa = self.kappa
          # set grey transmission function
          rad.Trans = rad.TransGrey
          # set flux function
          rad.LWFlux = rad.LWgrey

      elif rce_type == 'semigrey':
        print "doing semi-grey, 1 band RCE calculation"
        #semi-grey band data-- should be an arg, maybe
        semiGreyBandParams = Dummy()
        semiGreyBandParams.nu1 = 650.
        semiGreyBandParams.nu2 = 700.
        semiGreyBandParams.kappa = self.kappa
        rad.bandData = [semiGreyBandParams]
        # set trans function
        rad.Trans = rad.TransSemiGrey
        # set flux function
        rad.LWFlux = rad.LWtotal

      if self.data is None:
        self.data = dict()

      list_of_keys = ('p', 'T','Tad', 'flux', 'heat', 'q')
      for Tg in Tglist:

        # check if Tg is already in self.data
        # if it is, start from previous calculation
        if Tg in self.data:
          print 'found an old state for Tg = ' +  str(Tg)
          T = self.data[Tg]['T']
          self.Tad = self.data[Tg]['Tad']
          self.p = self.data[Tg]['p']
          self.q = self.data[Tg]['q']
          self.FRESH = False
        else:  # do a new one
          print 'no old state found for Tg = ' +  str(Tg) + ', starting a new one'
          self.FRESH = True
          #Use the following for the dry adiabat
          #Tad = Tg*(p/ps)**phys.air.Rcp #Change gas if desired

          #Use the following for the moist adiabat
          m = phys.MoistAdiabat(phys.water,phys.air)
          ps = self.ps + self.psat(Tg) #Total surface pressure for this temperature
          pl, Tad, molarCon, massCon = m(ps,Tg,p)
          T = 250*np.ones(len(p),np.Float) #Initialize constant T
          self.Tad = Tad
          self.p = p
          #------------Set the greenhouse gas concentration------------------------
          self.q = self.rh*massCon #water-like Oobleck?
          #q[:] = 10.e-4 #Makes optical depth 10 for Oobleck
          #------------------------------------------------------------------------
        # do some steps
        if self.FRESH:
          T,flux,heat = self.steps(T, Tg, 50,self.dtime)
          T,flux,heat = self.steps(T, Tg, 80,self.dtime, strat_do = True)
          T,flux,heat = self.steps(T, Tg, 200,self.dtime/10, strat_do = True)
          #T,flux,heat = self.steps(T,Tg, 100,self.dtime/50, strat_do = True)
        else:
          T,flux,heat = self.steps(T,Tg, 200,self.dtime/50, strat_do = True)
          T,flux,heat = self.steps(T,Tg, 200,self.dtime/100, strat_do = True)

        profiles = dict()
        for key, val in zip(list_of_keys, (self.p, T, self.Tad, flux, heat, self.q)):
            profiles[key] = val
        self.data[Tg] = profiles
      return

  def steps(self, T, Tg, nSteps, dtime, water_do = True, strat_do = False):
      '''step toward rce'''
      for i in range(nSteps):
          #Do smoothing
          if i%5 == 0 & i>0:
              for j in range(1,len(T)-1):
                  T[j] = .25*T[j-1] + .5*T[j] + .25*T[j+1]
                  self.q[j] = .25*self.q[j-1] + .5*self.q[j] + .25*self.q[j+1]
          #
          flux, heat = rad.LWFlux(self.p,T,Tg,self.q)
          dT = heat*dtime
          #Limit the temperature change per step
          dT = np.where(dT>5.,5.,dT)
          dT = np.where(dT<-5.,-5.,dT)
          #Midpoint method time stepping
          flux, heat = rad.LWFlux(self.p,T+.5*dT,Tg,self.q) # LW flux and heat for all bands
          dT = heat*dtime
          #Limit the temperature change per step
          dT = np.where(dT>5.,5.,dT)
          dT = np.where(dT<-5.,-5.,dT)
          T += dT
          #
          dTmax = max(abs(dT)) #To keep track of convergence
          #Uncomment next line to do hard convective adjustment
          T = np.where(T<self.Tad,self.Tad,T)
          #
          # update q following temperature changes
          if water_do:
            for j, ptot in enumerate(self.p):
              self.q[j] = self.rh*self.psat(T[j])/ptot*self.massrat
            if strat_do:
              stratT = self.find_tropopause(Tg, T = T, Tad = self.Tad)
              stratq = self.q[np.where(T==stratT)]
              #stratq = self.q[np.where(T>self.Tad+2)][-1] # mass mixing ratio at tropopause
              for j, ptot in enumerate(self.p):
                if T[j] > self.Tad[j]+0.1:
                  self.q[j] = stratq # fix stratospheric mass mixing ratio
                else:
                  break

          if i%50 == 0:
              print i,T[0],T[n/2],T[-2],dTmax/dtime
      return T,flux,heat

  def plot_data(self,x,y, flipy = True):
    '''plot keys x and y'''
    for Tg in self.data:
      plt.plot(self.data[Tg][x], self.data[Tg][y], label = str(Tg))
      if flipy:
        plt.ylim(self.data[Tg][y][-1], self.data[Tg][y][0])
    plt.legend()
    return

  def find_tropopause(self, Tg, T = None, Tad = None):
    '''find where the atmosphere stops following a moist adiabat'''
    if T is None:
      T = self.data[Tg]['T']
      Tad = self.data[Tg]['Tad']

    strat = T[np.where(T > Tad)]
    for thing in strat[-1::-1]:
      if thing < 250:
          return thing

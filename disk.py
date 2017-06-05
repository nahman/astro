import numpy as np
import constants as ct
import matplotlib.pyplot as plt
from scipy.misc import derivative
from scipy.integrate import quad

class Disk:
    Z=2e-2 #gas metallicity
    std_temp=260 #at 1AU

    def __init__(self, mode: "'pow' for power-law or 'exp' for exponential; 'mix' to use mixed law", tau_fr, alpha, t_init: 'in Myr' =0, atmos='clean', dep_time: 'in Myr' =5, dep_factor=1):
        self.mode=mode
        self.tau_fr=tau_fr
        self.alpha=alpha
        self.dep_factor=dep_factor
        self.atmos=atmos

        #self.dep_time=dep_time
        if mode=='mix':
            self.dep_time= (10*ct.AU_to_cm)**2/(alpha*self.c_s(10)*self.H(10))/ct.Myr_to_s #Note that 10AU is about where xrays carve a gap in the disk.
        else:
            self.dep_time=dep_time
        self.t_init=t_init #NOT LOGSPACE
        self.t_switch = None #used to store time at which we switch from power to exp law, assuming self.mixed is true.

        


    #Functions that change based on disk type
    
    def Sigma_gas(self,a):                       #not time dependent. gives initial surface density.
        return 100*Disk.Sigma_mmsn(a)*self.dep_factor
    
    def factor(self,t):
 
        if self.mode=='exp':
            return np.exp(-t/self.dep_time)   
        
        pow_factor=(t/self.t_init)**(-5/4)
        
        if self.mode=='mix' and pow_factor<=1e-2:
            
            if self.t_switch==None:
                self.t_switch=t

            return (self.t_switch/self.t_init)**(-5/4)*np.exp(-(t-self.t_switch)/self.dep_time)
             
        return pow_factor 
    
    def factor_dot(self,t):
        exp_factor_dot = -1/self.dep_time*self.factor(t)
        if self.mode=='exp':
            return exp_factor_dot
        pow_factor=(t/self.t_init)**(-5/4)
        if self.mode=='mix' and pow_factor<=1e-2:
            return exp_factor_dot
        return -5/4*1/(self.t_init)*(t/self.t_init)**(-9/4)

    def Sigma_gas_t(self,a,t):
        return self.factor(t)*self.Sigma_gas(a)

    def rho_gas(self,a,t):                             #volumetric density of gas. input AU, output CGS
        return self.Sigma_gas_t(a,t)/Disk.H(a)

    def P_gas(self,a,t):                               #Pressure of gas. input AU, output CGS. Analytically, 
        return self.rho_gas(a,t)*Disk.c_s(a)**2             #this is about 59.86 (a/1AU)^{-13/4}

    def diff_P_gas(self,a,t):
        return self.factor(t)*-15*(a/10)**-2.5*Disk.c_s(a)*Disk.Omega(a)+self.Sigma_gas_t(a,t)*-1.75*(ct.k*260*ct.G*ct.Msun/(ct.mu*ct.m_H*ct.AU_to_cm**3))**.5*a**-2.75
    
    def disk_mass(self,start,end,t):                    #input in AU, output in Mearth. mass of disk in between two radii
        return 4*ct.pi*1e4*ct.AU_to_cm**2*self.factor(t)/ct.Mearth*((end/10)**.5-(start/10)**.5)

    #Functions that work for all types of disks
    @staticmethod
    def temp(a):                            
        return Disk.std_temp*(a**-.5)
    @staticmethod
    def c_s(a):                                  #sound speed. input AU, output CGS                   
        return (ct.k*Disk.temp(a)/(ct.mu*ct.m_H))**.5
    @staticmethod
    def Omega(a):                               #keplerian particle velocity at a, in AU
        return (ct.G*ct.Msun/(a*ct.AU_to_cm)**3)**.5 
    @staticmethod
    def H(a):
        return Disk.c_s(a)/Disk.Omega(a)        #disk scale height, input in AU, output in cgs

    #maybe i should throw in a self to make it easier to call in planet
    @staticmethod
    def Sigma_mmsn(a):                          #minimum mass Solar nebula. A surface density. input AU, output cgs
        return (a/10)**(-1.5)


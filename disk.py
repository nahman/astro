import numpy as np
import constants as ct
import matplotlib.pyplot as plt
from scipy.misc import derivative
from scipy.integrate import quad

class Disk:
    Z=2e-2 #gas metallicity
    std_temp=260 #at 1AU

    def __init__(
                self, 
                mode: "type of gas depletion to use. 'pow' for power-law or 'exp' for exponential; 'mix' starts as power and becomes exp", 
                tau_fr: 'dimensionless stopping time', 
                alpha: 'from shakura sunyaev prescription for viscosity', 
                Sigma_0: 'initial gas disk density in 100s of sigma_mmsn'=1,
                atmos: "'clean' for clean envelope, 'dusty' for an envelope with dust grains. affects atm accretion"='clean', 
                dep_time: 'e-folding time of gas depletion under the exponential law. gets auto-calculated in mixed law. in Myr' =5, 
                t_init: 'in Myr. acts as a sort of normalization to the power law.'=1e-3, 
                ):

        self.mode=mode
        self.tau_fr=tau_fr
        self.alpha=alpha
        self.Sigma_0=Sigma_0
        self.atmos=atmos

        if mode=='mix':
            #Note that 10AU is about where xrays carve a gap in the disk. Thus we use 10AU as the 'disk edge'
            self.dep_time=(10*ct.AU_to_cm)**2/(alpha*self.c_s(10)*self.H(10))/ct.Myr_to_s 
        else:
            self.dep_time=dep_time
        
        self.t_init=t_init 
        
        self.t_switch = None #used to store time at which we switch from power to exp law, assuming self.mixed is true.

    
    #Functions that change based on disk type
    
    #gives initial surface density.
    def Sigma_gas(self,a):                       
        return 100*Disk.Sigma_mmsn(a)*self.Sigma_0
    
    #time dependent portion of Sigma_gas_t; multiplies with initial surface density
    def factor(self,t):
        if self.mode=='exp':
            return np.exp(-t/self.dep_time)        
        pow_factor=(t/self.t_init)**(-5/4)      
        if self.mode=='mix' and pow_factor<=1e-2:           
            if self.t_switch==None:
                self.t_switch=t
            return (self.t_switch/self.t_init)**(-5/4)*np.exp(-(t-self.t_switch)/self.dep_time)           
        return pow_factor 
    
    #time derivative of factor
    def factor_dot(self,t):
        exp_factor_dot = -1/self.dep_time*self.factor(t)
        if self.mode=='exp':
            return exp_factor_dot
        pow_factor=(t/self.t_init)**(-5/4)
        if self.mode=='mix' and pow_factor<=1e-2:
            return exp_factor_dot
        return -5/4*1/(self.t_init)*(t/self.t_init)**(-9/4)

    #time-evolving surface density of disk
    def Sigma_gas_t(self,a,t):
        return self.factor(t)*self.Sigma_gas(a)

    #volumetric density of gas. input AU, output CGS
    def rho_gas(self,a,t):                             
        return self.Sigma_gas_t(a,t)/Disk.H(a)
    
    #Pressure of gas. input AU, output CGS. 
    def P_gas(self,a,t):                               
        return self.rho_gas(a,t)*Disk.c_s(a)**2             
    
    #derivate of pressure wrt a
    def diff_P_gas(self,a,t):
        return self.factor(t)*-15*(a/10)**-2.5*Disk.c_s(a)*Disk.Omega(a)+self.Sigma_gas_t(a,t)*-1.75*(ct.k*260*ct.G*ct.Msun/(ct.mu*ct.m_H*ct.AU_to_cm**3))**.5*a**-2.75
    
    #mass of disk in between two radii. input in AU, output in Mearth. 
    def disk_mass(self,start,end,t):                    
        return 4*ct.pi*1e4*ct.AU_to_cm**2*self.factor(t)/ct.Mearth*((end/10)**.5-(start/10)**.5)


    #Functions that work for all types of disks
    
    #disk temp. input AU, output K. normalized as std_temp @ 1AU
    @staticmethod
    def temp(a):                            
        return Disk.std_temp*(a**-.5)
    
    #sound speed. input AU, output CGS   
    @staticmethod
    def c_s(a):                                                 
        return (ct.k*Disk.temp(a)/(ct.mu*ct.m_H))**.5
    
    #angular velocity of particles at distance a AU
    @staticmethod
    def Omega(a):                               
        return (ct.G*ct.Msun/(a*ct.AU_to_cm)**3)**.5 
    
    #disk scale height, input in AU, output in cgs
    @staticmethod
    def H(a):
        return Disk.c_s(a)/Disk.Omega(a)        
    
    #minimum mass Solar nebula. A surface density. input AU, output cgs
    @staticmethod
    def Sigma_mmsn(a):                          
        return (a/10)**(-1.5)


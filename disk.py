import numpy as np
import constants as ct


class Disk:
    Z=2e-2 #gas metallicity
    std_temp=260 #at 1AU

    def __init__(self,mode,tau_fr,alpha,t_init=0,dep_time=5,dep_factor=1):
        self.mode=mode
        self.tau_fr=tau_fr
        self.alpha=alpha
        self.dep_factor=dep_factor
        if mode=='exp':
            self.dep_time=5
        self.t_init=t_init

    #Functions that change based on disk type
    
    def Sigma_gas(self,a):                       #not time dependent. gives initial surface density.
        return 100*Disk.Sigma_mmsn(a)*self.dep_factor
    
    def factor(self,t):
        if self.mode=='exp':
            return np.exp(-t/self.dep_time)    
        return (t/self.t_init)**(-5/4) #I should move t_init somewhere else
    
    def factor_dot(self,t):
        if self.mode=='exp':
            return -1/self.dep_time*self.factor(t)
        return -5/4*1/(self.t_init)*(t/self.t_init)**(-9/4)

    def Sigma_gas_t(self,a,t):
        return self.factor(t)*self.Sigma_gas(a)
    def rho_gas(self,a):                             #volumetric density of gas. input AU, output CGS
        return self.Sigma_gas(a)/Disk.H(a)

    def P_gas(self,a):                               #Pressure of gas. input AU, output CGS. Analytically, 
        return self.rho_gas(a)*Disk.c_s(a)**2             #this is about 59.86 (a/1AU)^{-13/4}


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
    @staticmethod
    def diff_P_gas(a):                          #derivative of P_gas wrt a. Input in AU.
        return -195*a**(-17/4)

    #maybe i should throw in a self to make it easier to call in planet
    @staticmethod
    def Sigma_mmsn(a):                          #minimum mass Solar nebula. A surface density. input AU, output cgs
        return (a/10)**(-1.5)

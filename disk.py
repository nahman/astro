import numpy as np
import constants as ct
import matplotlib.pyplot as plt
from scipy.misc import derivative
from scipy.integrate import quad

class Disk:
    rho_solid=5.0
    Z=2e-2 #gas metallicity
    std_temp=260 #at 1AU
    gamma=1

    def __init__(
                self,    
                pebble_size: 'radius of solid pebbles in cm', 
                alpha: 'from shakura sunyaev prescription for viscosity', 
                t_s: 'timescale parameter for new Sigma gas',
                Sigma_0: 'initial gas disk density in 100s of sigma_mmsn'=1,
                atmos: "'clean' for clean envelope, 'dusty' for an envelope with dust grains. affects atm accretion"='clean', 
                const_Sigma: 'set to True for Sigma=3000 everywhere'=False,
                drag_mode: 'manually enforce epstein or stokes drag'='auto',
                fixed_tau: 'manually fix tau'=None
                ):

        self.pebble_size=pebble_size
        self.alpha=alpha
        self.Sigma_0=Sigma_0
        self.atmos=atmos
        self.t_s=t_s
        self.const_Sigma=const_Sigma
        self.drag_mode=drag_mode
        self.fixed_tau=fixed_tau

        self.R1=self.R_1_func(t_s)
        self.nu1=self.nu(self.R1)
        self.C=self.normalize(10*ct.m_jup)

        if not self.check_stbl(self.R1,0,True):
            raise Exception('unstable disk')
        
    #volumetric density of gas. input AU, output CGS
    def rho_gas(self,a,t):                             
        return self.new_Sigma_gas_t(a,t)/Disk.H(a)
    
    #Pressure of gas. input AU, output CGS. 
    def P_gas(self,a,t):                               
        return self.rho_gas(a,t)*Disk.c_s(a)**2             
    
    #derivate of pressure wrt a, cgs
    def diff_P_gas(self,a,t):
        gamma=Disk.gamma
        T=self.T(t)
        r=self.r(a)
        R1=self.R1

        ret=  1/R1*self.Sigma_0*self.C/(ct.pi*3*self.nu1)*T**(-(5/2-gamma)/(2-gamma))*(-gamma*r**(-gamma-1)*np.exp(-r**(2-gamma)/T)+r**-gamma*np.exp(-r**(2-gamma)/T)*(2-gamma)*-r**(1-gamma)/T)*self.c_s(a)*self.Omega(a) \
        + self.new_Sigma_gas_t(a,t)*-1.75*(ct.k*260*ct.G*ct.Msun/(ct.mu*ct.m_H*ct.AU_to_cm**3))**.5*a**-2.75
    
        return ret/ct.AU_to_cm

    #mass of disk in between two radii. input in AU, output in Mearth. 
    def disk_mass(self,start,end,t,show_acc=False):   

        f = lambda a: self.new_Sigma_gas_t(a,t)*2*ct.pi*a*ct.AU_to_cm**2/ct.Mearth

        if show_acc:
            return quad(f,start,end)     
        return quad(f,start,end)[0] 
    
    #mass accretion rate of star in Mearth/Myr
    def disk_mass_dot(self,t):
        gamma = Disk.gamma
        return self.C*self.T(t)**(-(5/2-gamma)/(2-gamma))*ct.AU_to_cm**2/ct.Mearth

    #viscosity (not molecular). input AU, output in AU^2/Myr
    def nu(self, a,cgs=False): 
        cgs_nu = self.alpha*self.c_s(a)*self.H(a)
        if cgs:
            return cgs_nu
        return cgs_nu*ct.AU_to_cm**-2*ct.Myr_to_s
    
    #calculate r, the unitless radial scale factor used in Hartmann 98
    def r(self,a):
        return a/self.R1
    
    #calculate T, the unitless time factor used in Hartmann 98
    def T(self,t):
        return t/self.t_s+1
    
    #calculate R1 for a give t_s (and implicitly, alpha)
    def R_1_func(self,t_s):
        gamma = Disk.gamma
        return t_s*3*(2-gamma)**2*self.alpha*ct.k*260/(ct.mu*ct.m_H)*(ct.AU_to_cm**3/(ct.G*ct.Msun))**.5*ct.AU_to_cm**-2*ct.Myr_to_s
    
    #sigma_gas model as in Hartmann 98
    def new_Sigma_gas_t(self,a,t):
        if self.const_Sigma:
            return 3000
        gamma = Disk.gamma
        ret = self.Sigma_0*self.C/(3*ct.pi*self.nu1*self.r(a)**gamma) * self.T(t)**(-(5/2-gamma)/(2-gamma))*np.exp(-self.r(a)**(2-gamma)/self.T(t))
        return ret
    
    #function for checking if the Toomre Q param is larger than one, and thus if the disk is grav stable
    def check_stbl(self,a,t,_bool=False):
        Q=self.c_s(a)*self.Omega(a)/(ct.pi*ct.G*self.new_Sigma_gas_t(a,t))
        if not _bool:
            return Q
        if Q>=1:
            return True
        return False

    #analytic expression for initial disk mass. does not factor in sigma_0
    def init_disk_mass(self):
        return 2/(3*self.nu1)*self.C*self.R1**2*ct.AU_to_cm**2/ct.Mearth
    
    #calculates C such that the total mass of the disk is some disk_mass
    def normalize(self,disk_mass):
        return disk_mass*ct.Mearth/ct.AU_to_cm**2*3*self.nu1/2/self.R1**2

    #Sigma gas divided by Sigma mmen
    def new_factor(self,a,t):
        ret = self.new_Sigma_gas_t(a,t)/self.Sigma_mmen(a)
        return ret

    #headwind factor: multiply this with keplerian velocity to find headwind. Look in OK12 for more info.
    def eta(self,a,t):                                 
        return -1/2*(Disk.c_s(a)/Disk.v_kep(a))**2*(a*ct.AU_to_cm/self.P_gas(a,t))*self.diff_P_gas(a,t) 

    #headwind velocity, cm/s
    def v_hw(self,a,t):                                
        return self.eta(a,t)*self.v_kep(a) 

    #solid particle scale height
    def H_p(self,a,t):
        return Disk.H(a)*np.sqrt(self.alpha/(self.alpha+self.tau(a,t)))
    
    def lmfp(self,a,t):
        return 1/((self.rho_gas(a,t)/(ct.mu*ct.m_H))*1e-15)

    #for dynamic tau, but fixed pebble size. Quadratic drag unimplemented
    def tau(self,a,t):
        if self.fixed_tau !=None:
            return self.fixed_tau
        #select a drag regime
        rho_gas=self.rho_gas(a,t)
        n=rho_gas/(ct.mu*ct.m_H)
        sigma=1e-15
        lmfp=1/(n*sigma)

        c_s=Disk.c_s(a)
        nu_molec=c_s*lmfp

        tau_ep=(self.rho_solid/rho_gas) * (self.pebble_size/c_s)*Disk.Omega(a)
        tau_st=1/(6*ct.pi)*(self.rho_solid/rho_gas)*self.pebble_size**2/nu_molec*Disk.Omega(a)
        
        #epstein
        if self.drag_mode=='epstein' or (self.drag_mode=='auto' and self.pebble_size<=9/4*lmfp):
            #print('epstein')
            return tau_ep
        #stokes
        else:
            #print('stokes')
            return tau_st



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
    
    #Keplerian velocty, cm/s
    @staticmethod
    def v_kep(a):                               
        return (ct.G*ct.Msun/(a*ct.AU_to_cm))**.5

    #disk scale height, input in AU, output in cgs
    @staticmethod
    def H(a):
        return Disk.c_s(a)/Disk.Omega(a)        
    
    #minimum mass Solar nebula. A surface density. input AU, output cgs
    @staticmethod
    def Sigma_mmsn(a):                     
        return (a/10)**(-1.5)
    
    #minimum mass extrasolar nebula
    @staticmethod
    def Sigma_mmen(a):
        return 4*10**5 * (a/.1)**-1.6
    


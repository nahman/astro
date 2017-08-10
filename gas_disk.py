import numpy as np
import constants as ct
from scipy.misc import derivative
from scipy.integrate import quad



class gas_disk:
    rho_solid=5.0
    Z=2e-2 #gas metallicity
    std_temp=260 #at 1AU
    gamma=1

    def __init__(
                self,     
                alpha: 'from shakura sunyaev prescription for viscosity', 
                t_s: 'timescale parameter for new Sigma gas',
                atmos: "'clean' for clean envelope, 'dusty' for an envelope with dust grains. affects atm accretion"='clean', 
                const_Sigma: 'set to True for Sigma=3000 everywhere'=False,
                gap_type: 'dc or c'='dc'
                ):

        self.alpha=alpha
        self.atmos=atmos
        self.t_s=t_s
        self.const_Sigma=const_Sigma

        self.R1=self.R_1_func(t_s)
        self.nu1=self.nu(self.R1)
        self.C=self.normalize(10*ct.m_jup)

        if not self.check_stbl(self.R1,0,True) and not const_Sigma:
            raise Exception('unstable disk')

        #dc (discontinuous gap): gap depletes to .1 unperturbed upon gap opening
        #c (continuous gap): gap depth evolves according to eqn 28 in DC13
        self.gap_type=gap_type
        
    def __enter__(self):
        return self
    
    def __exit__(self, exctype, excval, exc_trace):
        if excval != None:
            print('error: '+str(excval))
        return True


    #volumetric density of gas. input AU, output CGS
    def rho_gas(self,a,t):                             
        return self.new_Sigma_gas_t(a,t)/self.H(a)
    
    #Pressure of gas. input AU, output CGS. 
    def P_gas(self,a,t):                               
        return self.rho_gas(a,t)*self.c_s(a)**2             
    
    #derivate of pressure wrt a, cgs
    def diff_P_gas(self,a,t):
        gamma=gas_disk.gamma
        T=self.T(t)
        r=self.r(a)
        R1=self.R1

        ret=  1/R1*self.C/(ct.pi*3*self.nu1)*T**(-(5/2-gamma)/(2-gamma))*(-gamma*r**(-gamma-1)*np.exp(-r**(2-gamma)/T)+r**-gamma*np.exp(-r**(2-gamma)/T)*(2-gamma)*-r**(1-gamma)/T)*self.c_s(a)*self.Omega(a) \
        + self.new_Sigma_gas_t(a,t)*-1.75*(ct.k*260*ct.G*ct.Msun/(ct.mu*ct.m_H*ct.AU_to_cm**3))**.5*a**-2.75
    
        return ret/ct.AU_to_cm

    #mass accretion rate of star in Msun/yr
    def disk_mass_dot(self,t):
        gamma = gas_disk.gamma
        return self.C*self.T(t)**(-(5/2-gamma)/(2-gamma))*ct.AU_to_cm**2/ct.Msun* 1e-6

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
        gamma = gas_disk.gamma
        return t_s*3*(2-gamma)**2*self.alpha*ct.k*260/(ct.mu*ct.m_H)*(ct.AU_to_cm**3/(ct.G*ct.Msun))**.5*ct.AU_to_cm**-2*ct.Myr_to_s
    
    #sigma_gas model as in Hartmann 98



    def new_Sigma_gas_t(self,a,t,mp=None,gap=False):

        if self.const_Sigma:
            return 3000
        gamma = gas_disk.gamma
        sigma0 = self.C/(3*ct.pi*self.nu1*self.r(a)**gamma) * self.T(t)**(-(5/2-gamma)/(2-gamma))*np.exp(-self.r(a)**(2-gamma)/self.T(t))

        if gap:
            return .1*sigma0

        return sigma0

    #function for checking if the Toomre Q param is larger than one, and thus if the disk is grav stable
    def check_stbl(self,a,t,_bool=False):
        Q=self.c_s(a)*self.Omega(a)/(ct.pi*ct.G*self.new_Sigma_gas_t(a,t))
        if not _bool:
            return Q
        if Q>=1:
            return True
        return False

    #analytic expression for initial disk mass. 
    def init_disk_mass(self):
        return 2/(3*self.nu1)*self.C*self.R1**2*ct.AU_to_cm**2/ct.Mearth
    
    #calculates C such that the total mass of the disk is some disk_mass
    def normalize(self,disk_mass):
        return disk_mass*ct.Mearth/ct.AU_to_cm**2*3*self.nu1/2/self.R1**2

    #headwind factor: multiply this with keplerian velocity to find headwind. Look in OK12 for more info.
    def eta(self,a,t):                                 
        return -1/2*(self.c_s(a)/self.v_kep(a))**2*(a*ct.AU_to_cm/self.P_gas(a,t))*self.diff_P_gas(a,t) 

    #headwind velocity, cm/s
    def v_hw(self,a,t):                                
        return self.eta(a,t)*self.v_kep(a) 
    
    def lmfp(self,a,t):
        return 1/((self.rho_gas(a,t)/(ct.mu*ct.m_H))*1e-15)


    #Functions that work for all types of disks
    
    #disk temp. input AU, output K. normalized as std_temp @ 1AU
    @staticmethod
    def temp(a):                            
        return gas_disk.std_temp*(a**-.5)
    
    #sound speed. input AU, output CGS   
    @staticmethod
    def c_s(a,T=None):            
        if T !=None:
            return (ct.k*T/(ct.mu*ct.m_H))**.5                                      
        return (ct.k*gas_disk.temp(a)/(ct.mu*ct.m_H))**.5
    
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
        return gas_disk.c_s(a)/gas_disk.Omega(a)        
    
    #minimum mass Solar nebula. A surface density. input AU, output cgs
    @staticmethod
    def Sigma_mmsn(a):                     
        return (a/10)**(-1.5)
    
    #minimum mass extrasolar nebula
    @staticmethod
    def Sigma_mmen(a):
        return 4*10**5 * (a/.1)**-1.6
    
    @staticmethod
    def mach(a):
        return gas_disk.v_kep(a)/gas_disk.c_s(a)



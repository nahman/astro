import numpy as np
import constants as ct
import disk as d
import ID
import json
from scipy.misc import derivative
import matplotlib.pyplot as plt

def cubic(a2,a1,a0):
  #gives the real solution to the cubic equation
  #(Wolfram)

  sqrt = np.sqrt
  Q = (3*a1-a2**2)/9
  R = (9*a2*a1 -27*a0 -2*a2**3)/54
  D = Q**3 +R**2

  #now the number become complex
  S = (R+sqrt(D+0j))**(1./3)
  T = (R-sqrt(D+0j))**(1./3)

  #only the first solution is real; return this one
  x1 = -a2/3 +(S+T)
  #x2 = -a2/3 -(S+T)/2 +0.5j*sqrt(3)*(S-T)
  #x3 = -a2/3 -(S+T)/2 -0.5j*sqrt(3)*(S-T)
  x1=x1.real
  #print(x1)
  return x1

     
class Planet:
    
    C_acc=1.0
    GCR_fudge=2.8

    def __init__(
                self,
                m_c0: 'inital core mass in Mearth',
                m_a0: 'initial atm mass in Mearth',
                a0: 'inital orbital distance in AU',
                disk: d.disk,   
                case: '1: use type 1 migration till gap opening, then use type 2. 2: use type 1 till gas giant formation, then stop migrating.' =1,
                m_iso_mode: "set this to a mass in mearth to hard set pebble isolation mass. otherwise, set to 'auto' for automatic calculation via DM13" ='auto',
                reg: "'set', 'hyp', or '3b.' These regimes of core growth are laid out in OK12. choose 'auto' for automatic selection"  ='auto',
                const_core: 'sets time deriv of core mass to 0'=False,
                const_atm: 'sets time deriv of atm mass to 0'=False,
                const_a: 'sets time deriv of orbital distance to 0'=False,    
                corrected=True
                ):

        self.a=a0                       #core orbital radius
        self.m_core=m_c0                #core mass
        self.m_atm=m_a0                 #atmospher mass
        self.t=0                        #initial planet age in myr
        self.disk=disk                  #disk object that the planet is in
         
        #basic run parameters
        self.reg=reg                    #core accretion regime. set to 'set', '3b','hyp' or 'auto' for automatic selection
        self.m_iso_mode = m_iso_mode    #can be used to manually set an upperbound for planet mass, instead of dynamically calculating it
        self.const_core=const_core      #keeps core mass constant during run
        self.const_atm=const_atm        #keeps atmosphere mass constant during run
        self.const_a=const_a            #keeps orbital radius constant during run
        
        self.case=case                  #case 0: type 1 migration only ; 
                                        #case 1: type 1 -> type 2 upon gap opening ; 
                                        #case 2: type 1 -> type 2 upon Jupiter formation.
        
        if m_iso_mode=='auto':
            self.m_iso=self.gap_mass() 
        else:
            self.m_iso=m_iso_mode
        
        #the 'ID' object of the planet, used for easy identification. check ID.py for more info
        self.id=ID.Planet_ID(           
                    disk.pebble_size,
                    disk.alpha,
                    m_iso_mode,
                    disk.t_s,
                    disk.atmos,
                    a0,
                    case)
        
        #initializing a bunch of timestamps

        self.t_stop=None                #will be set to time at which pebble accretion ends
        self.t_jup=None                 #will be set to time of jupiter formation.
        self.t_death=None               #will be set to time of death, if necessary
        self.t_gap=None                 #will be set to time at which gap opens. I think functionally the same as t_stop
        
        #misc params
        self.toggled=False              #used to ensure certain things are printed only once. doesn't affect anything significant.
        self.died=False                 #boolean. did the planet die (via infall into the sun?)
        self.type2slow = 1              #premultiplies the type 2 migration rate. can be fudged to see the effect of faster or slower type 2.
        
        self.corrected=corrected        #'corrected' here means if the 3D correction factor given in OK12 is being taken into account.
                                        #the 3D correction factor accounts for the fact that the solid disk is not flat, but puffy, due to turbulence.

        self._mc_dot=None
        self._ma_dot=None
        self._a_dot=None

        self.t_dif=None

    def __enter__(self):
        return self
    
    def __exit__(self, exctype, excval, exc_trace):
        if excval != None:
            print('error: '+str(excval))
        return True

    #this function runs at every step while the planet's system of odes is being solved.
    def redefine(self,mc,ma,a,t):      
        #update the planet's stats
        self.m_core=mc
        self.m_atm=ma
        self.a=a
        self.t=t
        
        self._mc_dot=self.m_core_dot(t)

        if self.t_stop!=None:
            self.t_dif=t-self.t_stop

        
        self._ma_dot=self.m_atm_dot(t)
        self._a_dot=self.a_dot(t)

        if self.m_iso_mode=='auto':
            self.m_iso=self.gap_mass()

        if self.check_iso() and self.t_gap==None:
            self.t_gap=t
        elif not self.check_iso and self.t_gap!=None:
            self.t_gap=None

        if self.t_jup==None and self.m_atm+self.m_core>=ct.m_jup:
            self.t_jup=t
            if self.t_gap==None:
                self.t_gap=t
        elif self.t_jup!=None and self.m_atm+self.m_core<ct.m_jup:
            self.t_jup==None
        
        return self
    
    def state(self):
        return [self.m_core,self.m_atm,self.a]

    '''
    some basic functions
    '''   
    
    #angular velocity of planet in rad/s
    def Omega(self):
        return (ct.G*ct.Msun/(self.a*ct.AU_to_cm)**3)**.5

    #period of planet in days.
    def period(self): 
        return 2*np.pi/self.Omega() * (1/60) * (1/60) * (1/24)

    #bondi radius in cm
    def r_bondi(self):
        return ct.G*(self.m_atm+self.m_core)*ct.Mearth/self.disk.c_s(self.a)**2

    #hill radius, in cm
    def r_hill(self,mass=None,a=None):
        if mass == None and a == None:
            return self.a*ct.AU_to_cm*((self.m_core+self.m_atm)*ct.Mearth/(3*ct.Msun))**(1/3)
        return a*ct.AU_to_cm*(mass*ct.Mearth/(3*ct.Msun))**(1/3)

    #hill velocity, cm/s
    def v_hill(self):
        return self.r_hill()*self.Omega()
    
    #core radius, cm
    def r_core(self):
        return ((3*self.m_core*ct.Mearth)/(4*np.pi*self.disk.rho_solid))**(1/3)
    
    #core size in Hill units. 
    def alpha_core_alt(self):                   
        return self.r_core()/self.r_hill()
    
    def alpha_core(self):
        return 1e-3*(5/self.a)*(self.disk.rho_solid/3)**(-1/3)

    #Keplerian velocty, cm/s
    def v_kep(self):                               
        return (ct.G*ct.Msun/(self.a*ct.AU_to_cm))**.5
    
    #headwind factor: multiply this with keplerian velocity to find headwind. Look in OK12 for more info.
    def eta(self,t):                                 
        a=self.a
        return -1/2*(self.disk.c_s(a)/self.v_kep())**2*(a*ct.AU_to_cm/self.disk.P_gas(a,t))*self.disk.diff_P_gas(a,t) 

    #headwind velocity, cm/s
    def v_hw(self,t):        
        #print(self.eta(t)*self.v_kep() )                        
        return self.eta(t)*self.v_kep() 

    #OK12's 'headwind parameter.' Input AU, grams. Unitless.
    def zeta_w(self,t):                           
        return self.v_hw(t)/self.v_hill()
    
    #checks if the pebble isolation mass has been reached. 
    def check_iso(self):  
        return self.m_core+self.m_atm>=self.m_iso 

    '''
    Core Growth Functions
    '''
    #From OK12
    def impactB(self,t):
        zeta_w = self.zeta_w(t)
        St=self.disk.tau(self.a,t)

        al = self.alpha_core()

        bset=0
        bhyp=0
        b3b=0

        St_crit = min(12/zeta_w**3, 2)
        if St<=St_crit:
            bset = cubic(2*zeta_w/3, 0, -8*St) *np.exp(-(St/St_crit)**0.65)
        elif St>=St_crit and St>=zeta_w:
            b3b = (2*0.85*al**0.5 +1/St) *np.exp(-(0.7*zeta_w/St)**5)
        else:#intermediate case (give all)
            va = zeta_w/(1+St**2) *np.sqrt(1+4*St**2) #approach velocity
            bhyp = al *np.sqrt(1+6/(al*va**2))
            b3b = (2*0.85*al**0.5 +1/St) *np.exp(-(0.7*zeta_w/St)**5)
            bset = cubic(2*zeta_w/3, 0, -8*St) *np.exp(-(St/St_crit)**0.65)
        return [bset, bhyp, b3b]

    #From OK12
    def P_cols(self,t):

        """
        returns [bsets, vsets]

        returns the analalytical collision radius and approach velocity 
        in dimensionless units
        """
        #analytical collision rate
        #Approach velocity

        zeta_w = self.zeta_w(t)
        St=self.disk.tau(self.a,t)

        [bset, bhyp, b3b] = self.impactB(t)
        v3b = 3.2
        vhyp = zeta_w/(1+St**2) *np.sqrt(1+4*St**2)
        vset = 3*bset/2 +zeta_w 

        return [[bset,bhyp,b3b],[vset,vhyp,v3b]]
    
    #From OK12
    def P_col(self,t):
        """
        Obtain the (dimensionless) collision rate in the 2D regime:
            
            Pcol = 2R_col va
            (collision rate per unit sigma_d)
        """
        [barr, varr] = self.P_cols(t)
        Pcols = 2*np.array(barr)*np.array(varr)
        #pick the largest Pcol
        id = np.argmax(Pcols)
        '''
        if id==0:
            print('set',t)
        elif id==1:
            print('hyp',t)
        else:
            print('3b',t)
        '''
        factor = 1
        
        if self.corrected:
            factor = barr[id]*self.r_hill()/self.disk.H_p(self.a,t)

        #print(Pcols[id],barr[id],t)
        Pcol_ret,bcol_ret = Pcols[id]*min(1,factor) , barr[id]

        return [Pcol_ret,bcol_ret]

    #time derivative of core mass, in Mearth/Myr
    def m_core_dot(self,t):
        mdot = Planet.C_acc*self.P_col(t)[0]*self.disk.sigma_solid(self.a,t)*self.r_hill()*self.v_hill()*ct.Myr_to_s/ct.Mearth  
        
        if self.const_core or self.check_iso() or self.m_core>self.disk.init_solid_mass/ct.Mearth or self.t_gap!=None:
            mdot=0

        if mdot==0 and self.t_stop==None:
            self.t_stop=t
        elif self.t_stop!=None:
            self.t_stop==None

        #m_dot = Planet.C_acc*self.P_col(t)[0]*self.disk.sigma_solid(self.a,t)*self.r_hill()*self.v_hill()*ct.Myr_to_s/ct.Mearth  
        return mdot

    #accretion luminosity, erg/s
    def l_acc(self,t,cm=None):
        if cm!=None:
            self.m_core=cm

        self._mc_dot=self.m_core_dot(t)
        return ct.G*self._mc_dot*ct.Mearth/ct.Myr_to_s*self.m_core*ct.Mearth/self.r_core()
   
    '''
    Atmosphere growth functions
    '''
    def rho_rcb(self,t,tstart=None):
        #GCR=self.m_atm/self.m_core
        rcore=self.r_core()
        
        if self.disk.atmos=='clean':
            g=4/3
            T=self.disk.temp(self.a)
            GCR=self.GCR_clean(t,tstart)
            grad=ct.grad_adb_clean
        else:
            g=1.2
            T=ct.T_rcb
            #GCR=self.GCR_dusty(t)
            grad=ct.grad_adb
        
        Rb_rcb = ct.G*self.m_core*ct.Mearth/self.disk.c_s(self.a,T=T)**2
        exp = 1/(g-1)
        num = GCR*self.m_core*ct.Mearth
        den = 4*ct.pi*(grad*Rb_rcb)**exp*rcore**(3-exp)
        #print(num/den)
        return num/den
        
    def opacity(self,t,tstart=None):
        if self.disk.atmos=='clean':
            T=self.disk.temp(self.a)
            ret = 1e-5*(self.rho_rcb(t,tstart)/1e-6)**.6*(T/100)**2.2*self.disk.Z/.02
        else:
            T=ct.T_rcb
            ret = .03*(self.rho_rcb(t,tstart)/1e-4)**.5*(T/2500)**7.5*self.disk.Z/.02
        return ret
   
    #cooling luminosity, erg/s
    def l_cool(self,t,mc=None,tstart=None):

        if mc!=None:
            self.m_core=mc
        T = None
        GCR = self.m_atm/self.m_core
        mcore = self.m_core*ct.Mearth
        opacity = self.opacity(t,tstart)
        
        rho_rcb = self.rho_rcb(t,tstart=tstart)
        if self.disk.atmos =='clean':
            gamma=4/3
            T = self.disk.temp(self.a)
            GCR = self.GCR_clean(t,tstart=tstart)
        else:
            gamma=1.2
            T = ct.T_rcb
            #GCR = self.GCR_dusty(t)
        
        num = 64*ct.pi*ct.G*(1+GCR)*mcore*ct.sb*T**3*ct.mu_rcb*ct.m_H*(gamma-1)/gamma
        den = 3*ct.k*rho_rcb*opacity

        if den == 0:
            return np.inf

        return num/den
        
    def GCR_clean(self,t): #input in Myr, Mearth
        a = self.a

        if self.t_gap!=None and self.t>=self.t_gap:
            _gap=True
        else:
            _gap=False

        sigma_dep = .37*(self.disk.new_Sigma_gas_t(a,t,gap=_gap)/(self.disk.Sigma_mmen(a)*.03))**.12*(372/self.disk.temp(a))**1.5
        F_clean = sigma_dep*(.02/d.disk.Z)**.4*(ct.grad_adb_clean/.25)**2.2*(ct.mu_rcb/2.37)**2.2
        
        return F_clean*(self.t_dif*10)**.4*(self.m_core/5)

    def GCR_clean_dot(self,t): 
        
        #avoid non-analyticity @t=0 errors
        if self.t_dif<=1e-6: 
            return 0

        ret=self.GCR_clean(t)*4*(10*self.t_dif)**-1  
        return ret
    
    def m_atm_dot_clean(self,t):
        
        m_core=self.m_core
        m_dot = self.GCR_clean_dot(t)*m_core 
        
        return m_dot
    
    def GCR_dusty(self,t):    
        t_dif = self.t_dif
        a = self.a
        
        if self.t_gap!=None and self.t>=self.t_gap:
            _gap=True
        else:
            _gap=False

        sigma_dep = .03 * (self.disk.new_Sigma_gas_t(a,t,gap=_gap)/(self.disk.Sigma_mmen(a)*.03))**.12
        F_dusty = sigma_dep*(2500/ct.T_rcb)**4.8*(.02/self.disk.Z)**.4*(ct.grad_adb/.17)**3.4*(ct.mu_rcb/2.37)**3.4 
        return F_dusty*(10*t_dif)**.4*(self.m_core/5)**1.7

    def GCR_dusty_dot(self,t): 
        t_dif=self.t_dif
        if t_dif<=1e-6:
            return 0
        ret = self.GCR_dusty(t)*4*(10*t_dif)**-1
        return ret

    def m_atm_dot_dusty(self,t):

        m_core=self.m_core
        m_dot = self.GCR_dusty_dot(t)*m_core

        return m_dot
    
    def m_atm_dot(self,t):
        
        m_core=self.m_core
        m_atm=self.m_atm

        if self.const_atm or m_atm+m_core>=ct.m_jup or self.t_stop==None:
            return 0

        #Once the GCR hits 50%, runaway accretion kicks in. 
        if m_atm/m_core>=.5 and self.disk.disk_mass(t,self.a)>=ct.m_jup: 
            return 5e6 #this number needs to be more variable. a fixed high number seems to make odeint more prone to errors

        mdot=None
        if self.disk.atmos=='dusty':
            mdot=self.m_atm_dot_dusty(t)
        elif self.disk.atmos=='clean':
            mdot=self.m_atm_dot_clean(t)
        else:
            raise Exception('check your clean/dusty input')

        if mdot<=1e-10:
            return 0
        return mdot
 
    '''
    Migration functions
    '''

    #type 2 migration rate
    def a_dot_visc(self,t):
        loc_disk_m = (self.a*ct.AU_to_cm)**2*self.disk.new_Sigma_gas_t(self.a,t)
        #print(loc_disk_m/ct.Mearth/ct.m_jup)
        m_p = (self.m_core+self.m_atm) * ct.Mearth
        return -self.disk.alpha*(ct.k/(ct.mu*ct.m_H)*self.disk.std_temp*(ct.AU_to_cm)**.5/(ct.G*ct.Msun)**.5)*ct.Myr_to_s/ct.AU_to_cm*min(1,loc_disk_m/m_p)*self.type2slow

    #normalization torque
    def Gamma_0(self,t): 
        a=self.a
        _mp = self.m_core+self.m_atm
        ret = (_mp*ct.Mearth/ct.Msun)**2*(self.disk.H(a)/(a*ct.AU_to_cm))**-2*self.disk.new_Sigma_gas_t(a,t)*(a*ct.AU_to_cm)**4*self.disk.Omega(a)**2
        return ret
    
    #total torque
    def Gamma_tot(self,t):
        return -(1.36+.62*ct.beta_sigma+.43*ct.beta_t)*self.Gamma_0(t)

    #time derivative of planet's orbital distance. inputs in AU, Mearth, Myr. 
    def a_dot(self,t):
        
        if self.const_a:
            return 0 

        m_core=self.m_core*ct.Mearth
        m_atm=self.m_atm*ct.Mearth
        cgs_a_dot = 2*self.Gamma_tot(t)/((m_core+m_atm)*(ct.G*ct.Msun)**.5*(self.a*ct.AU_to_cm)**-.5)
        ret = cgs_a_dot*ct.Myr_to_s/ct.AU_to_cm
        
        #case 1
        if self.case==1:
            if (self.t_gap!=None and t>=self.t_gap) or self.check_iso(): 
                #for a smoother transition, uncomment the two lines below. Not the higher the sharpness factor, the quicker the transition occurs.
                #sharpness=10000
                return self.a_dot_visc(t)#*(1-np.exp(-sharpness*(t-self.t_gap)))+ret*np.exp(-sharpness*(t-self.t_gap))
            
        #case 2
        elif self.case==2:
            if self.t_jup!=None and t>=self.t_jup:
                return 0

        return ret
        
    #from Duffel/McFayden 13
    def gap_mass(self):
        a=self.a
        disk_mach = self.disk.mach(a)
        m_sh=(ct.Msun/ct.Mearth)/(disk_mach)**3
        return m_sh*17*np.sqrt(self.disk.alpha*disk_mach)

    #same as above but input params are explicitly set.
    def gap_mass_alt(self,a):

        disk_mach = self.disk.mach(a)
        m_sh=(ct.Msun/ct.Mearth)/(disk_mach)**3
        return m_sh*17*np.sqrt(self.disk.alpha*disk_mach)


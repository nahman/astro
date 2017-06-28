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
                disk: d.Disk,   
                case: '1: use type 1 migration till gap opening, then use type 2. 2: use type 1 till gas giant formation, then stop migrating.' =1,
                m_iso_mode: "set this to a mass in mearth to hard set pebble isolation mass. otherwise, set to 'auto' for automatic calculation via DM13" ='auto',
                reg: "'set', 'hyp', or '3b.' These regimes of core growth are laid out in OK12. choose 'auto' for automatic selection"  ='auto',
                const_core: 'sets time deriv of core mass to 0'=False,
                const_atm: 'sets time deriv of atm mass to 0'=False,
                const_a: 'sets time deriv of orbital distance to 0'=False,    
                corrected=False
                ):

        self.a=a0                       #core orbital radius
        self.m_core=m_c0                #core mass
        self.m_atm=m_a0                 #atmospher mass
        self.t=1e-3                     #initial planet age in myr
        self.disk=disk                  #disk object that the planet is in
        self.m_iso_mode = m_iso_mode    

        if m_iso_mode=='auto':
            #kinda risky calling a function of this sort inside __init__. needs an alternative solution
            self.m_iso=self.gap_mass() 
        else:
            self.m_iso=m_iso_mode

        self.reg=reg
        self.const_core=const_core
        self.const_atm=const_atm
        self.const_a=const_a
        self.case=case
        
        self.tstop=0
                         #will be set to the time at which the core stops growing
        self.astop=None                 #will be set to orbital distance at which core stops growing
        

        self.id=ID.Planet_ID(           #the 'ID' object of the planet, used for easy identification. check ID.py for more info
                    disk.pebble_size,
                    disk.alpha,
                    m_iso_mode,
                    disk.t_s,
                    disk.atmos,
                    case)
        
        self.t_jup=None                 #will be set to time of jupiter formation.
        self.toggled=False              #used to ensure certain things are printed only once. doesn't affect anything significant.
        self.died=False
        self.t_death=None

        self.corrected=corrected

    def redefine(self,mc,ma,a,t):       #this function runs at every step while the planet's system of odes is being solved.

        #update the planet's stats
        self.m_core=mc
        self.m_atm=ma
        self.a=a
        self.t=t

        if self.m_iso_mode=='auto':
            self.m_iso=self.gap_mass()

        return self
    
    def state(self):
        return [self.m_core,self.m_atm,self.a]

    def reset(self,y0,t0):
        self.m_core=y0[0]
        self.m_atm=y0[1]
        self.a=y0[2]
        self.t=t0
        return self


    '''
    some basic functions
    '''   

    def Omega(self):
        return (ct.G*ct.Msun/(self.a*ct.AU_to_cm)**3)**.5

    #period of planet in days.
    def period(self): 
        return 2*np.pi/self.Omega() * (1/60) * (1/60) * (1/24)

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

    #bondi radius, cm. c_s should be in cgs
    def r_bondi(self,c_s: 'in cm/s'):                           
        return ct.G*(self.m_atm+self.m_core)*ct.Mearth/(c_s**2)

    #Keplerian velocty, cm/s
    def v_kep(self):                               
        return (ct.G*ct.Msun/(self.a*ct.AU_to_cm))**.5
    
    #headwind factor: multiply this with keplerian velocity to find headwind. Look in OK12 for more info.
    def eta(self,t):                                 
        a=self.a
        return -1/2*(d.Disk.c_s(a)/self.v_kep())**2*(a*ct.AU_to_cm/self.disk.P_gas(a,t))*self.disk.diff_P_gas(a,t) 

    #headwind velocity, cm/s
    def v_hw(self,t):                                
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

        Pcol_ret,bcol_ret = Pcols[id]*min(1,factor) , barr[id]

        return [Pcol_ret,bcol_ret]

    #time derivative of core mass, in Mearth/Myr
    def m_core_dot(self,t):

        if self.const_core:
            return 0
        
        if self.check_iso():
            if not self.toggled:
                if np.around(self.m_core,5)==np.around(self.m_iso,5) or self.tstop==0:
                    self.tstop=t
                if self.astop==None:
                    self.astop=self.a
                self.toggled=True
            return 0
          
        m_dot = Planet.C_acc*self.P_col(t)[0]*self.disk.Sigma_mmsn(self.a)*self.r_hill()*self.v_hill()*ct.Myr_to_s/ct.Mearth  
  
        return m_dot

    '''
    Atmosphere growth functions
    '''
    def GCR_clean(self,t): #input in Myr, Mearth
        t_dif=t-self.tstop
        
        a = self.a
        F_clean = .1*(Planet.GCR_fudge/2.8)*(200/self.disk.temp(a))**1.5*(.02/d.Disk.Z)**.4*(ct.grad_adb_clean/.25)**2.2*(ct.mu_rcb/2.37)**2.2
        
        return F_clean*(1000*t_dif)**.4*(self.m_core/5)*self.disk.new_factor(a,t)**.12

    def GCR_clean_dot(self,t): 
        
        t_dif=t-self.tstop #avoid non-analyticity @t=0 errors
        if t_dif<=0:
            return 0

        ret=self.GCR_clean(t)*400*(1000*t_dif)**-1  

        return ret
    
    def m_atm_dot_clean(self,t):
        m_core=self.m_core

        #I shouldn't need to calculate the latter half since it will be 0.
        m_dot = self.GCR_clean_dot(t)*m_core #+self.GCR_clean(t)*self.m_core_dot(t)
        
        return m_dot
    
    def GCR_dusty(self,t):    
        t_dif = t-self.tstop

        F_dusty = .16*(Planet.GCR_fudge/3)*(2500/ct.T_rcb)**4.8*(.02/self.disk.Z)**.4*(ct.grad_adb/.17)**3.4*(ct.mu_rcb/2.37)**3.4 #factors from eqn 20 of FC14, plus an efolidng scaling on nebular density
        GCR=F_dusty*t_dif**.4*(self.m_core/5)**1.7*self.disk.new_factor(self.a,t)**.12 
        return GCR

    def GCR_dusty_dot(self,t): 
        t_dif = t-self.tstop
        ret = self.GCR_dusty(t)*.4*t_dif**-1

    def m_atm_dot_dusty(self,t):

        m_core=self.m_core
        m_atm=self.m_atm

        m_dot = self.GCR_dusty_dot(t)*m_core+self.GCR_dusty(t)*self.m_core_dot(t)

        return m_dot
    
    def m_atm_dot(self,t):
        
        m_core=self.m_core
        m_atm=self.m_atm

        if self.const_atm or not self.check_iso() or self.disk.new_Sigma_gas_t(self.a,t)<=1e-7 or m_atm>=ct.m_jup:
            return 0

        #Once the GCR hits 50%, runaway accretion kicks in. 
        if m_atm/m_core>=.5 and self.disk.disk_mass(self.a,100,t)>=ct.m_jup: 
            return 1e5 

        mdot=None
        if self.disk.atmos=='dusty':
            mdot=self.m_atm_dot_dusty(t)
        elif self.disk.atmos=='clean':
            mdot=self.m_atm_dot_clean(t)
        else:
            raise Exception('check your clean/dusty input')
        return mdot
 
    '''
    Migration functions
    '''
    def a_dot_visc(self,t):
        loc_disk_m = (self.a*ct.AU_to_cm)**2*self.disk.new_Sigma_gas_t(self.a,t)
        m_p = (self.m_core+self.m_atm) * ct.Mearth
        return -self.disk.alpha*(ct.k/(ct.mu*ct.m_H)*260*(ct.AU_to_cm)**.5/(ct.G*ct.Msun)**.5)*ct.Myr_to_s/ct.AU_to_cm * min(1,loc_disk_m/m_p)

    #normalization torque
    def Gamma_0(self,t): 
        a=self.a
        m = self.m_core+self.m_atm
        ret = (m*ct.Mearth/ct.Msun)**2*(self.disk.H(a)/(a*ct.AU_to_cm))**-2*self.disk.new_Sigma_gas_t(a,t)*(a*ct.AU_to_cm)**4*self.disk.Omega(a)**2
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

        #case 1
        if self.case==1:
            if self.check_iso(): 
                #print(self.m_iso)
                #return 0 #turning of viscous migration
                return self.a_dot_visc(t)
            
        #case 2
        elif self.case==2:
            if self.m_atm+self.m_core>=ct.m_jup:
                return 0

        cgs_a_dot= 2*self.Gamma_tot(t)/((m_core+m_atm)*(ct.G*ct.Msun)**.5*(self.a*ct.AU_to_cm)**-.5)
        return cgs_a_dot*ct.Myr_to_s/ct.AU_to_cm
        
    #from Fung, Shi, Chiang 2014
    def check_gap(self):
        a=self.a
        m=(self.m_atm+self.m_core)*ct.Mearth
        q = m/ct.Msun               #planet to Sun mass ratio
        h_r= self.disk.H(a)/(a*ct.AU_to_cm)   #disk aspect ratio

        if 1e-4<=q<5e-3:
            if .14*(q/1e-3)**-2.16*(ct.alpha_ss/1e-2)*(h_r/.05)**6.61 <= self.gap_param:
                return True

            return False

        elif 5e-3<=q<=1e-2:
            if (q/5e-3)**-1.00*(ct.alpha_ss/1e-2)**1.26*(h_r/.05)**5.12 <= self.gap_param:
                return True

            return False

        elif q<1e-4:
            return False            #I assume low mass planets don't form gaps
    
        else:
            #if the planet has gotten this big, something has gone wrong.
            raise Exception('check_gap failed. the mass of the planet exceeds upper bounds set for the calculation of gap formation')
        
    #from Duffel/McFayden 13
    def gap_mass(self):
        a=self.a
        disk_mach = a*ct.AU_to_cm/self.disk.H(a)
        m_sh=(ct.Msun/ct.Mearth)/(disk_mach)**3
        return m_sh*17*np.sqrt(self.disk.alpha*disk_mach)


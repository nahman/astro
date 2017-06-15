import numpy as np
import constants as ct
import disk as d
import ID
import json
from scipy.misc import derivative

class Planet:
    rho_solid=5.0
    C_acc=1.0
    GCR_fudge=2.0
    Pcol_fudge=3.0
    smooth=1e10
    gap_param=1.0

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
                ):

        self.a=a0                       #core orbital radius
        self.m_core=m_c0                #core mass
        self.m_atm=m_a0                 #atmospher mass
        self.t=disk.t_init              #planet age
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
        
        self.tstop=None                 #will be set to the time at which the core stops growing
        self.astop=None                 #will be set to orbital distance at which core stops growing
        

        self.id=ID.Planet_ID(           #the 'ID' object of the planet, used for easy identification. check ID.py for more info
                    disk.tau_fr,
                    disk.alpha,
                    m_iso_mode,
                    disk.t_s,
                    disk.atmos,
                    case)
        
        self.t_jup=None                 #will be set to time of jupiter formation.
        self.toggled=False              #used to ensure certain things are printed only once. doesn't affect anything significant.


    def redefine(self,mc,ma,a,t):       #this function runs at every step while the planet's system of odes is being solved.

        #update the planet's stats
        self.m_core=mc
        self.m_atm=ma
        self.a=a
        self.t=t
        #print('a: '+str(a))
        #print('Sigma: '+str(self.disk.new_Sigma_gas_t(a,t)))
        #print('m_c: '+str(mc))
        #print('m_a: '+str(ma))
        #recalculate the new pebble isolation mass (which is dependent on a), if necessary
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
    def r_hill(self):
        return self.a*ct.AU_to_cm*((self.m_core+self.m_atm)*ct.Mearth/(3*ct.Msun))**(1/3)
    
    #hill velocity, cm/s
    def v_hill(self):
        return self.r_hill()*self.Omega()
    
    #core radius, cm
    def r_core(self):
        return ((3*self.m_core*ct.Mearth)/(4*np.pi*self.rho_solid))**(1/3)
    
    #core size in Hill units. 
    def alpha_core(self):                   
        return self.r_core()/self.r_hill()

    #bondi radius, cm. c_s should be in cgs
    def r_bondi(self,c_s: 'in cm/s'):                           
        return ct.G*(self.m_atm+self.m_core)*ct.Mearth/(c_s**2)

    #Keplerian velocty, cm/s
    def v_kep(self):                               
        return (ct.G*ct.Msun/(self.a*ct.AU_to_cm))**.5
    
    #headwind factor: multiply this with keplerian velocity to find headwind. Look in OK12 for more info.
    def eta(self,t):                                 
        a=self.a
        return -1/2*(d.Disk.c_s(a)/(self.v_kep()**2))**2*(a*ct.AU_to_cm/self.disk.P_gas(a,t))*self.disk.diff_P_gas(a,t) 

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
    def Pcol(self,t):
        zeta = self.zeta_w(t)
        r_H = self.r_hill()
        v_H = self.v_hill()
        alpha = self.alpha_core()
        tau_fr=self.disk.tau_fr

        #settling calc
        coeff = [1,2/3*zeta,0,-8*tau_fr]
        coeff=np.nan_to_num(coeff)
        roots = np.roots(coeff)
        b_set_arr=roots.real[abs(roots.imag)<1e-5]
        b_set=max(b_set_arr)                                #get real positive root
        
        v_a_set=3/2*b_set+zeta
        f_set = np.exp(-(tau_fr/min(12/zeta**3,2))**.65)    #smoothing factor
        P_col_set = 2*b_set*f_set*v_a_set 

        #hyperbolic calc
        v_a_hyp=zeta*(1+4*tau_fr**2)**.5/(1+tau_fr)**2
        b_hyp=alpha*(1+6/(alpha*v_a_hyp**2))**.5
        P_col_hyp = 2*b_hyp*v_a_hyp

        #three body calc
        b_3b=1.7*alpha+1/tau_fr
        v_a_3b=3.2
        f_3b=np.exp(-(.7*zeta/tau_fr)**5)                   #smoothing factor


        P_col_3b = 2*b_3b*f_3b*v_a_3b

        #check if we should force a regime or let the max pick a regime for us.

        if self.reg=='auto':
            ret = max(P_col_set,P_col_hyp,P_col_3b)

            #I just often use this for checking
            '''
            if ret==P_col_set:
                print('set: '+str(P_col_set)+'3b: '+str(P_col_3b))
                
            elif ret==P_col_hyp:
                print('hyp')
            else:
                print('3b: '+str(P_col_3b)+'set: '+str(P_col_set))
            '''

            return ret
        elif self.reg =='set':
            return P_col_set
        elif self.reg =='hyp':
            return P_col_hyp
        elif self.reg=='3b':
            return P_col_3b
        else:
            raise Exception('error in Pcol. no regime with such name')

    #time derivative of core mass, in Mearth/Myr
    def m_core_dot(self,t):

        if self.const_core:
            return 0
        
        if self.check_iso():
            if self.tstop==None:
                self.tstop=t
            if self.astop==None:
                self.astop=self.a
            return 0
          
        m_dot = Planet.C_acc*self.Pcol(t)/self.Pcol_fudge*self.disk.Sigma_mmsn(self.a)*self.r_hill()*self.v_hill()*ct.Myr_to_s/ct.Mearth  
  
        return m_dot

    '''
    Atmosphere growth equations
    '''
    def GCR_clean(self,t): #input in Myr, Mearth
        a = self.a
        F_clean = .1*(Planet.GCR_fudge/2.8)*(200/self.disk.temp(a))**1.5*(.02/d.Disk.Z)**.4*(ct.grad_adb_clean/.25)**2.2*(ct.mu_rcb/2.37)**2.2

        return F_clean*(1000*t)**.4*(self.m_core/5)*self.disk.new_factor(a,t)**.12

    def GCR_clean_dot(self,t): 
        ret=self.GCR_clean(t)*(400*(1000*t)**-1+(self.m_core/5)**-1*(self.m_core_dot(t)/5)+.12*self.disk.new_factor(self.a,t)**-1*self.disk.new_factor_dot(self.a,t))
        return ret
    
    def m_atm_dot_clean(self,t):
        m_core=self.m_core
        m_atm=self.m_atm

        #I shouldn't need to calculate the latter half since it will be 0.
        m_dot = self.GCR_clean_dot(t)*m_core#+self.GCR_clean(t)*self.m_core_dot(t)
        
        return m_dot
    
    def GCR_dusty(self,t):              
        F_dusty = .16*(Planet.GCR_fudge/3)*(2500/ct.T_rcb)**4.8*(.02/self.disk.Z)**.4*(ct.grad_adb/.17)**3.4*(ct.mu_rcb/2.37)**3.4 #factors from eqn 20 of FC14, plus an efolidng scaling on nebular density
        GCR=F_dusty*t**.4*(self.m_core/5)**1.7*self.disk.new_factor(self.a,t)**.12 
        return GCR

    def GCR_dusty_alt(self,m_core,t):
        F_dusty = .16*(Planet.GCR_fudge/3)*(2500/ct.T_rcb)**4.8*(.02/self.disk.Z)**.4*(ct.grad_adb/.17)**3.4*(ct.mu_rcb/2.37)**3.4 #factors from eqn 20 of FC14, plus an efolidng scaling on nebular density
        GCR=F_dusty*t**.4*(m_core/5)**1.7*self.disk.new_factor(self.a,t)**.12 
        return GCR

    def GCR_dusty_dot(self,t): 
        return self.GCR_dusty(t)*((.4*t**-1)+1.7/5*(self.m_core/5)**-1*self.m_core_dot(t)+.12*self.disk.new_factor(self.a,t)**-1*self.disk.new_factor_dot(self.a,t))

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
            return 1e8 #*factor #this expression is really large, then goes to 0 when m_atm hits m_jup.
        
        mdot=None
        if self.disk.atmos=='dusty':
            mdot=self.m_atm_dot_dusty(t)
        else:
            mdot=self.m_atm_dot_clean(t)
        
        #print(mdot)
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
                #return 0 #turning of viscous migration
                return self.a_dot_visc(t)
            
        #case 2
        elif self.case==2:
            if self.m_atm+self.m_core>=ct.m_jup:
                return 0

        cgs_a_dot= 2*self.Gamma_tot(t)/((m_core+m_atm)*(ct.G*ct.Msun)**.5*(self.a*ct.AU_to_cm)**-.5)
        return cgs_a_dot*ct.Myr_to_s/ct.AU_to_cm
        
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

    def gap_mass(self):
        a=self.a
        disk_mach = a*ct.AU_to_cm/self.disk.H(a)
        m_sh=(ct.Msun/ct.Mearth)/(disk_mach)**3
        return m_sh*17*np.sqrt(self.disk.alpha*disk_mach)
        

    def check_gap_alt(self):
        return self.m_atm+self.m_core>=self.gap_mass()
    


import numpy as np
import constants as ct
import disk as d
import ID
import json
from collections import namedtuple


class Planet:
    rho_solid=5.0
    C_acc=1.0
    GCR_fudge=2.0
    Pcol_fudge=3.0
    core_smooth=1e+10
    atm_smooth=1e+10

    def __init__(self,m_c0,m_a0,a0,m_iso,disk,reg='auto',const_core=False,const_atm=False,const_a=False,use_iso=True):
        self.a=a0
        self.m_core=m_c0
        self.m_atm=m_a0
        self.m_iso=m_iso
        self.reg=reg
        self.const_core=const_core
        self.const_atm=const_atm
        self.const_a=const_a
        self.tstop=None
        self.use_iso=use_iso
        self.id=ID.Planet_ID(disk.tau_fr,disk.alpha,m_iso,disk.mode,disk.t_init)
        self.disk=disk  #is composition a good idea?

    def redefine(self,mc,ma,a):
        self.m_core=mc
        self.m_atm=ma
        self.a=a
    
    def state(self):
        return [self.m_core,self.m_atm,self.a]

    '''
    some basic functions
    '''   
    def Omega(self):
        return (ct.G*ct.Msun/(self.a*ct.AU_to_cm)**3)**.5

    def r_hill(self):
        return self.a*ct.AU_to_cm*((self.m_core+self.m_atm)*ct.Mearth/(3*ct.Msun))**(1/3)
    
    def v_hill(self):
        return self.r_hill()*self.Omega()
        
    def r_core(self):
        return ((3*self.m_core)/(4*np.pi*self.rho_solid))**(1/3)
    
    def alpha_core(self):                   #core size in Hill units. input AU, grams.
        return self.r_core()/self.r_hill()

    def r_bondi(self,c_s):                           #bondi radius. input mearths, AU, output CGS
        return ct.G*(self.m_atm+self.m_core)*ct.Mearth/(c_s**2)

    def check_iso(self):  
        if self.m_core+self.m_atm>=self.m_iso: #Supposed to be both, right?
            return True
        return False

    def v_kep(self):                               #Keplerian velocty. input in AU, output in cgs
        return (ct.G*ct.Msun/(self.a*ct.AU_to_cm))**.5
    
    def eta(self):                                 #headwind factor: multiply this with keplerian velocity to find headwind. Input in AU. Unitless.
        a=self.a
        return -1/2*(self.disk.c_s(a)/(self.v_kep()**2))**2*(a*ct.AU_to_cm/self.disk.P_gas(a))*self.disk.diff_P_gas(a) 

    def v_hw(self):                                #headwind velocity. input AU, output cgs.
        return self.eta()*self.v_kep() 

    def zeta_w(self):                            #OK12's 'headwind parameter.' Input AU, grams. Unitless.
        return self.v_hw()/self.v_hill()

    '''
    Core Growth Functions
    '''
    def Pcol(self):
        zeta = self.zeta_w()
        r_H = self.r_hill()
        v_H = self.v_hill()
        alpha = self.alpha_core()
        tau_fr=self.disk.tau_fr

        #settling calc
        coeff = [1,2/3*zeta,0,-8*tau_fr]
        coeff=np.nan_to_num(coeff)
        roots = np.roots(coeff)
        b_set=roots.real[abs(roots.imag)<1e-5][0]
        v_a_set=3/2*b_set+zeta
        f_set = np.exp(-(tau_fr/min(12/zeta**3,2))**.65) #smoothing factor

        P_col_set = 2*b_set*f_set*v_a_set 

        #hyperbolic calc
        v_a_hyp=zeta*(1+4*tau_fr**2)**.5/(1+tau_fr)**2
        b_hyp=alpha*(1+6/(alpha*v_a_hyp**2))**.5

        P_col_hyp = 2*b_hyp*v_a_hyp

        #three body calc
        b_3b=1.7*alpha+1/tau_fr
        v_a_3b=3.2
        f_3b=np.exp(-(.7*zeta/tau_fr)**5) #smoothing factor
        
        P_col_3b = 2*b_3b*f_3b*v_a_3b

        #check if we should force a regime or let the max do its thing.

        if self.reg=='auto':
            return max(P_col_set,P_col_hyp,P_col_3b)
        elif self.reg =='set':
            return P_col_set
        elif self.reg =='hyp':
            return P_col_hyp
        elif self.reg=='3b':
            return P_col_3b
        else:
            raise Exception('error in Pcol. no regime with such name')

    def m_core_dot(self,t):
        if self.const_core:
            return 0

        m_dot = Planet.C_acc*self.Pcol()/self.Pcol_fudge*self.disk.Sigma_mmsn(self.a)*self.r_hill()*self.v_hill()*ct.Myr_to_s/ct.Mearth 

        if self.check_iso():
            if self.tstop==None:
                self.tstop=t
            return m_dot*np.exp(-Planet.core_smooth*np.abs(self.m_core-self.m_iso))
        return m_dot

    '''
    Atmosphere growth equations
    '''

    def GCR_dusty(self,t):              
        F_dusty = .16*(Planet.GCR_fudge/3)*(2500/ct.T_rcb)**4.8*(.02/self.disk.Z)**.4*(ct.grad_adb/.17)**3.4*(ct.mu_rcb/2.37)**3.4 #factors from eqn 20 of FC14, plus an efolidng scaling on nebular density
        GCR=F_dusty*t**.4*(self.m_core/5)**1.7*self.disk.factor(t)**.12 
        return GCR

    def GCR_dusty_dot(self,t): 
        return self.GCR_dusty(t)*((.4*t**-1)+1.7/5*(self.m_core/5)**-1*self.m_core_dot(t)+.12*self.disk.factor(t)*-1*self.disk.factor_dot(t))

    def m_atm_dot_dusty(self,t):
        m_core=self.m_core
        m_atm=self.m_atm

        if self.const_atm or m_core<=1:
            return 0
        
        m_dot = self.GCR_dusty_dot(t)*m_core+self.GCR_dusty(t)*self.m_core_dot(t)

        #Once the GCR hits 50%, runaway accretion kicks in. 
        if m_atm/m_core>=.5: 
            factor=1
            if m_atm>=ct.m_jup:
                factor=0
            
            return m_dot*1e8*factor+(1-factor)*np.exp(-1e6*np.abs(m_atm-ct.m_jup)) #USE EXPONENTIALS TO DRIVE THE MASS TO MJUP
            
            #note that the exponential is being used as a smoothing factor to keep derivatives continuous, 

        return m_dot
    
    '''
    Migration functions
    '''
    def a_dot_visc(self):
        return-self.disk.alpha*(ct.k/(ct.mu*ct.m_H)*260*(ct.AU_to_cm)**.5/(ct.G*ct.Msun)**.5)*ct.Myr_to_s/ct.AU_to_cm

    #normalization torque 
    def Gamma_0(self,t): 
        a=self.a
        m = self.m_core+self.m_atm
        return (m*ct.Mearth/ct.Msun)**2*(self.disk.H(a)/(a*ct.AU_to_cm))**-2*self.disk.Sigma_gas_t(a,t)*(a*ct.AU_to_cm)**4*self.disk.Omega(a)**2

    #total torque
    def Gamma_tot(self,t):
        return -(1.36+.62*ct.beta_sigma+.43*ct.beta_t)*self.Gamma_0(t)

    #time derivative of planet's orbital distance. inputs in AU, Mearth, Myr. 
    def a_dot(self,t): 
        if self.const_a:
            return 0 

        m_core=self.m_core*ct.Mearth
        m_atm=self.m_atm*ct.Mearth
        
        if self.use_iso:
            if self.check_iso(): 
                return 0 #turning of viscous migration
                #return self.a_dot_visc()

            cgs_a_dot= 2*self.Gamma_tot(t)/((m_core+m_atm)*(ct.G*ct.Msun)**.5*(self.a*ct.AU_to_cm)**-.5)

            
            return cgs_a_dot*ct.Myr_to_s/ct.AU_to_cm
        
        raise Exception('Only iso mass implemented right now')

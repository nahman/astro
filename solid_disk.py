import numpy as np
import disk as d    
import constants as ct
import scipy.interpolate as interp
from scipy.integrate import quad
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class solid_disk(d.Disk):

    def __init__(self,    
                pebble_size: 'radius of solid pebbles in cm', 
                alpha: 'from shakura sunyaev prescription for viscosity', 
                t_s: 'timescale parameter for new Sigma gas',
                sigma_edge_0: 'initial surface density at edge of solid disk',
                Sigma_0: 'initial gas disk density in 100s of sigma_mmsn'=1,
                atmos: "'clean' for clean envelope, 'dusty' for an envelope with dust grains. affects atm accretion"='clean', 
                const_Sigma: 'set to True for Sigma=3000 everywhere'=False,
                drag_mode: 'manually enforce epstein or stokes drag'='auto',
                fixed_tau: 'manually fix tau'=None

                ):
        d.Disk.__init__(self,pebble_size,alpha,t_s,Sigma_0,atmos,const_Sigma,drag_mode,fixed_tau)
        self.sigma_edge_0=sigma_edge_0

        def diff_system(a,t):
            if a<=1e-3:
                a=1e-3
                return 0
            a_dot=-self.v_drift(a,t)
            return a_dot
        
        t=np.logspace(-3,1,1000)
        sol = odeint(diff_system,self.R1,t)[:,0]
        self.edge_data=(t,sol)

    #input in AU, Myr. you can choose units of output
    def v_drift(self,a,t,cgs=False):
        if a <=1e-3:
            return 0
        tau=self.tau(a,t)
        cgs_v=2*self.v_hw(a,t)*tau/(1+tau**2)
        if cgs:
            return cgs_v
        return cgs_v * ct.Myr_to_s/ct.AU_to_cm
    
    #input in Myr output in AU
    def edge_pos(self,t):
        f = interp.interp1d(self.edge_data[0],self.edge_data[1],kind='linear',fill_value=0)
        return f(t)
    
    #input in Myr, output in cgs (or whatever units sigma_edge_0 is in)
    def sigma_edge(self,t):
        r_edge=self.edge_pos(t)
        if r_edge<=1e-3:
            return 0

        return self.sigma_edge_0*self.R1/r_edge
    
    #input in Myr, output in cgs
    def m_dot(self,t):
        r_edge=self.edge_pos(t)
        if r_edge<=1e-3:
            return 0
        return self.sigma_edge(t)*2*ct.pi*r_edge*ct.AU_to_cm*self.v_drift(r_edge,t,True)
    
    def sigma_solid(self,a,t):
        if a >=self.edge_pos(t) or a<=1e-3:
            return 0
        ret=self.m_dot(t)/(2*ct.pi*a*ct.AU_to_cm*self.v_drift(a,t,True))

        return ret



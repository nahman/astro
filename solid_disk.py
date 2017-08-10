import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
from scipy.integrate import odeint, quad
from mpl_toolkits.mplot3d import Axes3D

import constants as ct
import gas_disk as gd


class disk(d.Disk):

    def __init__(self,    
                pebble_size: 'radius of solid pebbles in cm', 
                alpha: 'from shakura sunyaev prescription for viscosity', 
                t_s: 'timescale parameter for new Sigma gas',
                disk_edge: "how big the disk is: set to 'R1' for gas disk R1",
                sigma_1: 'initial surface density at 1AU of solid disk',
                Sigma_0: 'initial gas disk density in 100s of sigma_mmsn'=1,
                atmos: "'clean' for clean envelope, 'dusty' for an envelope with dust grains. affects atm accretion"='clean', 
                const_Sigma: 'set to True for Sigma=3000 everywhere'=False,
                drag_mode: 'manually enforce epstein or stokes drag'='auto',
                fixed_tau: 'manually fix tau'=None
                ):
        gd.gas_disk.__init__(self,pebble_size,alpha,t_s,Sigma_0,atmos,const_Sigma,drag_mode)

        self.pebble_size = pebble_size
        self.fixed_tau = fixed_tau

        if disk_edge=='R1':
            self.disk_edge = self.R1
        else:
            self.disk_edge = disk_edge
        
        self.sigma_1 = sigma_1

        def diff_system(a,t):
            if a<=1e-2:
                a=1e-2
                return 0
            a_dot=-self.v_drift(a,t)
            return a_dot
        
        t=np.logspace(-3,1,1000)
        sol = odeint(diff_system,self.disk_edge,t)[:,0]
        self.edge_data=(t,sol)

        f = lambda a: (a)**-1*sigma_1 * 2 * ct.pi * a * ct.AU_to_cm**2
        self.init_solid_mass = quad(f,1e-2,self.disk_edge)[0]
    
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


    #input in AU, Myr. you can choose units of output
    def v_drift(self,a,t,cgs=False):
        if a <=1e-2:
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
    
    #input in Myr, output in cgs (or whatever units sigma_edge is in)
    def sigma_solid(self,a,t):
        edge = self.edge_pos(t)
        if a>edge:
            return 0
        a=float(a)
        f = lambda a: a**-1 * 2 * ct.pi * a * ct.AU_to_cm**2
        C = quad(f,1e-2,edge)[0]
        
        return a**-1*self.init_solid_mass/C 
    
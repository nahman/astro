import numpy as np
import scipy.interpolate as interp
from scipy.integrate import odeint, quad

import constants as ct
import gas_disk as gd


#this disk inherits all methods from gas_disk, and introduces a time evolving solid component
class disk(gd.gas_disk):

    def __init__(self,    
                pebble_size: 'radius of solid pebbles in cm', 
                alpha: 'from shakura sunyaev prescription for viscosity', 
                t_s: 'timescale parameter for new Sigma gas',
                disk_edge: "how big the disk is: set to 'R1' for gas disk R1" = 'R1',
                #sigma_1: 'initial surface density at 1AU of solid disk',
                atmos: "'clean' for clean envelope, 'dusty' for an envelope with dust grains. affects atm accretion"='clean', 
                const_Sigma: 'set to True for Sigma=3000 everywhere'=False,
                drag_mode: 'manually enforce epstein or stokes drag'='auto',
                fixed_tau: 'manually fix tau'=None,
                gap_type: 'dc or c'='dc',
                s_mass: 'mass of solids in disk, mearth'=30
                ):
        gd.gas_disk.__init__(self,alpha,t_s,atmos,const_Sigma,gap_type)

        self.pebble_size = pebble_size
        self.fixed_tau = fixed_tau
        self.drag_mode=drag_mode

        if disk_edge=='R1':
            self.disk_edge = self.R1
        else:
            self.disk_edge = disk_edge

        def diff_system(a,t):
            if a<=1e-2:
                a=1e-2
                return 0
            a_dot=-self.v_drift(a,t)
            return a_dot
        
        t=np.logspace(-6,1,1000)
        sol = odeint(diff_system,self.disk_edge,t)[:,0]      

        a_exp = np.floor(np.log10(self.disk_edge))
        accuracy=5
        i=0

        while min(sol)<0 or max(sol)>self.disk_edge*1.5:
        
            if i>=5:
                raise Exception('unable to solve disk edge equation')
            
            new_edge = np.random.randint(-99,99)*(10**(a_exp-accuracy))+self.disk_edge
            print('disk wiggling: '+str(new_edge))
            sol = odeint(diff_system,new_edge,t)[:,0] 
            _min=min(sol)
            _max=max(sol)
            i+=1
                        
        self.edge_data=(t,sol)

        self.init_solid_mass = s_mass*ct.Mearth 
        self.t_edge=None
    
    #for dynamic tau, but fixed pebble size. Quadratic drag unimplemented
    def tau(self,a,t):
        if self.fixed_tau !=None:
            return self.fixed_tau

        #select a drag regime
        rho_gas=self.rho_gas(a,t)
        n=rho_gas/(ct.mu*ct.m_H)
        sigma=1e-15
        lmfp=1/(n*sigma)

        c_s=self.c_s(a)
        nu_molec=c_s*lmfp

        tau_ep=(self.rho_solid/rho_gas) * (self.pebble_size/c_s)*self.Omega(a)
        tau_st=4/9*(self.rho_solid/rho_gas)*self.pebble_size**2/nu_molec*self.Omega(a)
        
        #epstein
        if self.drag_mode=='epstein' or (self.drag_mode=='auto' and self.pebble_size<=9/4*lmfp):
            #print('epstein',t)
            return tau_ep
        #stokes
        else:
            #print('stokes',t)
            return tau_st
    
    #solid particle scale height
    def H_p(self,a,t):
        return self.H(a)*np.sqrt(self.alpha/(self.alpha+self.tau(a,t)))

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
        f = interp.interp1d(self.edge_data[0],self.edge_data[1],kind='linear',fill_value=1e-2,bounds_error=False)
        return f(t)
    
    #input in Myr, output in cgs (or whatever units sigma_edge is in)
    def sigma_solid(self,a,t):
        edge = self.edge_pos(t)

        if np.abs(a-edge)<=1e-6:
            return 0

        if a>edge:
            if self.t_edge==None:
                self.t_edge=t
            return 0

        a=float(a)
        f = lambda a: a**-1 * 2 * ct.pi * a * ct.AU_to_cm**2
        C = quad(f,1e-2,edge)[0]
        
        return a**-1*self.init_solid_mass/C 

    #mass of disk in between two radii. input in AU, output in Mearth. 
    def disk_mass(self,t,start,end='auto',show_acc=False,mode:"'gas' or 'solid'"='gas'):   

        if mode=='gas':
            sigma=self.new_Sigma_gas_t
            if end=='auto':
                end=self.R1*2
        elif mode=='solid':
            sigma=self.sigma_solid
            if end=='auto':
                end=self.edge_pos(t) 
            if end<=start:
                return 0
        else:
            raise Exception('check disk mass mode')


        #use unperturbed sigma_gas
        f = lambda a: sigma(a,t)*2*ct.pi*a*ct.AU_to_cm**2/ct.Mearth

        if show_acc:
            return quad(f,start,end)     
        
        ret,acc = quad(f,start,end)

        return ret

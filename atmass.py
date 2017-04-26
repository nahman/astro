import constants as ct
import coremass
import basic_equations as eqn 
import numpy as np

'''
Assumptions:
-gas accretion begins once the core hits 1 earth mass
-at GCR .5, runaway accretion kicks in, making the planet jupiter-massed. Right now I just set the accretion rate to 0.
-e folding sigma_gas
'''

#expressions related to a clean atmosphere. need to be reworked.
'''
def GCR_clean(a,m_core,t): #input in Myr, Mearth
    F_clean = .1*(ct.fudge/2.8)*(200/eqn.temp(a))**1.5*(.02/ct.Z)**.4*(ct.grad_adb/.25)**2.2*(ct.mu_rcb/2.37)**2.2
    return F_clean*(1000*t)**.4*(m_core/5)*eqn.e_fold_factor(t)**.12

def GCR_clean_dot(a,m_core,m_core_dot,t): 
    return GCR_clean(a,m_core,t)*(400*t**-1+(m_core/5)**-1*(m_core_dot/5)+.12*eqn.e_fold_factor(t)**-1*eqn.e_fold_factor_dot(t))

def m_atm_dot_clean(a,m_core,m_core_dot,m_atm,t):
    
    if m_atm/m_core>.5: 
        return 0

    return GCR_clean_dot(a,m_core,m_core_dot,t)*m_core+GCR_clean(a,m_core,t)*m_core_dot 

'''
#Gas-Core mass ratio. input in Myr, Mearth, AU
def GCR_dusty(a,m_core,t):              
    F_dusty = .16*(ct.fudge/3)*(2500/ct.T_rcb)**4.8*(.02/ct.Z)**.4*(ct.grad_adb/.17)**3.4*(ct.mu_rcb/2.37)**3.4 #factors from eqn 20 of FC14, plus an efolidng scaling on nebular density
    GCR=F_dusty*t**.4*(m_core/5)**1.7*eqn.Sigma_factor(t)**.12 
    return GCR

#Time derivative of GCR_dusty
def GCR_dusty_dot(a,m_core,m_core_dot,t): 
    return GCR_dusty(a,m_core,t)*((.4*t**-1)+1.7/5*(m_core/5)**-1*m_core_dot+.12*eqn.Sigma_factor(t)*-1*eqn.Sigma_factor_dot(t))

#time derivative of atmospheric mass, derived using GCR_dusty (assuming dusty atmosphere). inputs are in AU, Mearth, Mearth/Myr, Mearth, Myr.
def m_atm_dot_dusty(a,m_core,m_core_dot,m_atm,t):
    if ct.const_atm or m_core<=1:
        return 0
    
    m_dot = GCR_dusty_dot(a,m_core,m_core_dot,t)*m_core+GCR_dusty(a,m_core,t)*m_core_dot

    #Once the GCR hits 50%, runaway accretion kicks in. Right now, I just force this function to 0. 
    if m_atm/m_core>=.5: 
        
        return m_dot*np.exp(-ct.atm_smooth*(m_atm/m_core-.5)) #this is a place holder. need to make m_planet m_jup.
        
        #note that the exponential is being used as a smoothing factor to keep derivatives continuous, 
        #but also to force it to 0 quickly. odeint likes continuity.

    return m_dot




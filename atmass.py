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

def GCR_dusty(a,m_core,t): #input in Myr, Mearth, AU
    F_dusty = .16*(ct.fudge/3)*(2500/ct.T_rcb)**4.8*(.02/ct.Z)**.4*(ct.grad_adb/.17)**3.4*(ct.mu_rcb/2.37)**3.4 #factors from eqn 20 of FC14, plus an efolidng scaling on nebular density
    GCR=F_dusty*t**.4*(m_core/5)**1.7*eqn.e_fold_factor(t)**.12 
    return GCR

def GCR_dusty_dot(a,m_core,m_core_dot,t): 
    return GCR_dusty(a,m_core,t)*((.4*t**-1)+1.7/5*(m_core/5)**-1*m_core_dot+.12*eqn.e_fold_factor(t)*-1*eqn.e_fold_factor_dot(t))

def m_atm_dot_dusty(a,m_core,m_core_dot,m_atm,t):
    if ct.const_atm or m_core<=1:
        return 0
    
    if m_atm/m_core>=.5: 
        return (GCR_dusty_dot(a,m_core,m_core_dot,t)*m_core+GCR_dusty(a,m_core,t)*m_core_dot)*np.exp(-1000*(m_atm/m_core-.5)) #this is a place holder. I should figure out a way to signal m_p to be set to m_jup, or I should break the calculation here and solve another system of eqns where m_p=m-jup
        #note that the exponential is being used as a smoothign factor to keep my derivatives continuous
        
    return GCR_dusty_dot(a,m_core,m_core_dot,t)*m_core+GCR_dusty(a,m_core,t)*m_core_dot




import constants as ct
import numpy as np
from scipy.misc import derivative



#functions for getting various astro parameters

def r_core(m_core): #input, output cgs
    return ((3*m_core)/(4*ct.pi*ct.rho_solid))**(1/3)

def alpha_core(a,m_core): #input AU, grams, output hill units
    return r_core(m_core)/r_hill(a,m_core)

def temp(a): #AU, K
    return 260*(a**-.5)

def c_s(a): #AU, CGS
    return (ct.k*temp(a)/(ct.mu*ct.m_H))**.5

def r_bondi(m,a): #grams, AU
    return ct.G*m/(c_s(a)**2)

def Sigma_mmsn(a): #input in AU
    return (a/10)**(-1.5)

def Sigma_gas(a): #AU
    return 100*Sigma_mmsn(a)

def Sigma_gas_efold(a,t): #input in Myrs
    return Sigma_gas(a)*np.exp(-t/ct.depletion_time)

def Sigma_gas_efold_dot(a,t): #time derivative of sigma gas efold
    return -1/ct.depletion_time*Sigma_gas_efold(a,t)

def Omega(a): #input in AU. Keplerian. output rad/s
    return (ct.G*ct.Msun/(a*ct.AU_to_cm)**3)**.5

def H(a): #disk scale height, inputs in AU
    return c_s(a)/Omega(a)

def rho_gas(a): #input AU, output CGS
    return Sigma_gas(a)/H(a)

def P_gas(a): #analytically, this is abotu 59.86 (a/1AU)^{-13/4}
    return rho_gas(a)*c_s(a)**2

def diff_P_gas(a): #derivative of P_gas wrt a. Input in AU.
    return -195*a**(-17/4)

def r_hill(a,m): #input in au, gram. output in cm.
    return a*ct.AU_to_cm*(m/(3*ct.Msun))**(1/3)

def v_hill(a,m):
    return r_hill(a,m)*Omega(a)

def v_kep(a):
    return (ct.G*ct.Msun/(a*ct.AU_to_cm))**.5

def eta(a): #headwind factor: multiply this with keplerian velocity to find headwind
    return -1/2*(c_s(a)/(v_kep(a)**2))**2*(a*ct.AU_to_cm/P_gas(a))*diff_P_gas(a) 

def v_hw(a): #headwind velocity
    return eta(a)*v_kep(a) 

def zeta_w(a,m): #OK12's 'headwind parameter'
    return v_hw(a)/v_hill(a,m)


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
    return 100*Sigma_mmsn(a)*ct.depletion_factor

def e_fold_factor(t): #input in Myrs
    if ct.const_gas:
        return 1
    return np.exp(-t/ct.depletion_time)

def Sigma_gas_efold(a,t): #input in Myrs
    return Sigma_gas(a)*e_fold_factor(t)

def e_fold_factor_dot(t):
    if ct.const_gas:
        return 0
    return -1/ct.depletion_time*e_fold_factor(t)

def Sigma_gas_efold_dot(a,t): #time derivative of sigma gas efold
    return Sigma_gas(a)*e_fold_factor_dot(a,t)

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


#Condition checking for gap opening (from Eugene's paper)

def check_gap(a,m):
    q = m/ct.Msun

    h_r= H(a)/(a*ct.AU_to_cm) #disk aspect ratio

    if 1e-4<=q<5e-3:

        if .14*(q/1e-3)**-2.16*(ct.alpha_ss/1e-2)*(h_r/.05)**6.61 <= ct.gap_param:
            return True
        return False

    elif 5e-3<=q<=1e-2:

        if (q/5e-3)**-1.00*(ct.alpha_ss/1e-2)**1.26*(h_r/.05)**5.12 <= ct.gap_param:
            return True
        return False
    elif q<1e-4:
        return False #I assume low mass planets don't form gaps
    else:

        raise Exception('check_gap failed. the mass of the planet exceeds upper bounds set for the calculation of gap formation')

def bin_check_gap(a,m):
    if check_gap(a,m):
        return 0
    return 1
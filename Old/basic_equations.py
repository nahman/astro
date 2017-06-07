import constants as ct
import numpy as np
from scipy.misc import derivative


'''
functions for getting various astro parameters
'''

def r_core(m_core):                         #core radius, assuming some density. input, output cgs
    return ((3*m_core)/(4*ct.pi*ct.rho_solid))**(1/3)

def alpha_core(a,m_core):                   #core size in Hill units. input AU, grams.
    return r_core(m_core)/r_hill(a,m_core)

def temp(a):                                #temp of disk. takes AU, K
    return 260*(a**-.5)

def c_s(a):                                 #sound speed. input AU, output CGS
    return (ct.k*temp(a)/(ct.mu*ct.m_H))**.5

def r_bondi(m,a):                           #bondi radius. input grams, AU, outbut CGS
    return ct.G*m/(c_s(a)**2)

def Omega(a):                               #Keplerian angular velocity. input in AU, output rad/s
    return (ct.G*ct.Msun/(a*ct.AU_to_cm)**3)**.5

def H(a):                                   #disk scale height, input in AU, output in cgs
    return c_s(a)/Omega(a)

def rho_gas(a):                             #volumetric density of gas. input AU, output CGS
    return Sigma_gas(a)/H(a)

def P_gas(a):                               #Pressure of gas. input AU, output CGS. Analytically, 
    return rho_gas(a)*c_s(a)**2             #this is about 59.86 (a/1AU)^{-13/4}

def diff_P_gas(a):                          #derivative of P_gas wrt a. Input in AU.
    return -195*a**(-17/4)

def r_hill(a,m):                            #Hill radius. input in AU, gram. output in cm.
    return a*ct.AU_to_cm*(m/(3*ct.Msun))**(1/3)

def v_hill(a,m):
    return r_hill(a,m)*Omega(a)             #Hill velocity. input in AU, gram, output in cgs

def v_kep(a):                               #Keplerian velocty. input in AU, output in cgs
    return (ct.G*ct.Msun/(a*ct.AU_to_cm))**.5

def eta(a):                                 #headwind factor: multiply this with keplerian velocity to find headwind. Input in AU. Unitless.
    return -1/2*(c_s(a)/(v_kep(a)**2))**2*(a*ct.AU_to_cm/P_gas(a))*diff_P_gas(a) 

def v_hw(a):                                #headwind velocity. input AU, output cgs.
    return eta(a)*v_kep(a) 

def zeta_w(a,m):                            #OK12's 'headwind parameter.' Input AU, grams. Unitless.
    return v_hw(a)/v_hill(a,m)


'''
Important Gas Parameters
'''

def Sigma_mmsn(a):                          #minimum mass Solar nebula. A surface density. input AU, output cgs
    return (a/10)**(-1.5)

def Sigma_gas(a):                           #gas surface density of disk. input AU, output cgs
    return 100*Sigma_mmsn(a)

def power_factor(t):                        #power scaling factor that is multiplied with sigma_gas to give gas 
    return (t/10**ct.t_init)**(-5/4)        #surface density the appropriate scaling with time. input Myr.

def Sigma_gas_power(a,t):                     #Sigma_gas, but evolves with time according to a power rule. input Myr.
    return Sigma_gas(a)*power_factor(t)

def power_factor_dot(t):                    #time derivative of power_factor.
    return -5/4*1/(10**ct.t_init)*(t/10**ct.t_init)**(-9/4)

def Sigma_gas_power_dot(a,t):               #time derivative of sigma_gas_power
    return Sigma_gas(a)*power_factor_dot(t)

def e_fold_factor(t):                       #alternative factor that is multiplied with sigma_gas to give gas
    return np.exp(-t/ct.depletion_time)     #surface density the appropriate scaling with time. Efolding instead of power law. input Myr.

def Sigma_gas_efold(a,t):                   #Sigma_gas, but evolves with time according to e-fold timescale. input Myr.
    return Sigma_gas(a)*e_fold_factor(t)

def e_fold_factor_dot(t):                   #time derivative of e_fold_factor
    return -1/ct.depletion_time*e_fold_factor(t)

def Sigma_gas_efold_dot(a,t):               #time derivative of sigma_gas_efold
    return Sigma_gas(a)*e_fold_factor_dot(t)

def Sigma_factor(t):                         #unified sigma-multiplying factor that makes decision based on constants.py
    if ct.const_gas:
        return 1
    if ct.Sigma_pow:
        return power_factor(t)
    return e_fold_factor(t)

def Sigma_factor_dot(t):                    #time derivative of sigma_factor
    if ct.const_gas:
        return 0
    if ct.Sigma_pow:
        return power_factor_dot(t)
    return e_fold_factor_dot(t)

def Sigma_gas_t(a,t):                       #unified expression for time-evolving sigma_gas. 
    return Sigma_gas(a)*Sigma_factor(t)


'''
Gap checking functions
'''

#This function comes from FC14, eqns 12 and 14

def check_gap(a,m):
    q = m/ct.Msun               #planet to Sun mass ratio
    h_r= H(a)/(a*ct.AU_to_cm)   #disk aspect ratio

    if 1e-4<=q<5e-3:
        if .14*(q/1e-3)**-2.16*(ct.alpha_ss/1e-2)*(h_r/.05)**6.61 <= ct.gap_param:
            return True

        return False

    elif 5e-3<=q<=1e-2:
        if (q/5e-3)**-1.00*(ct.alpha_ss/1e-2)**1.26*(h_r/.05)**5.12 <= ct.gap_param:
            return True

        return False

    elif q<1e-4:
        return False            #I assume low mass planets don't form gaps
   
    else:
        #if the planet has gotten this big, something has gone wrong.
        raise Exception('check_gap failed. the mass of the planet exceeds upper bounds set for the calculation of gap formation')

'''
simple checking for isolation mass
'''

def check_iso(m):  #gram input
    if m/ct.Mearth>=ct.iso_mass:
        return True
    return False

print(Sigma_gas(1)*power_factor(10))

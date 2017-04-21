'''
This stuff is way too fast. The planet falls into the sun before a gap can form. It's causing expressions in coremass, atmass to basically go to infinity.
'''

import constants as ct
import basic_equations as eqn 
from coremass import check_gap

def Gamma_0(a,m_core,m_atm,t): #normalization torque (inputs in AU, GRAMS, myr)
    m = m_core+m_atm
    return (m/ct.Msun)**2*(eqn.H(a)/(a*ct.AU_to_cm))**-2*eqn.Sigma_gas_efold(a,t)*(a*ct.AU_to_cm)**4*eqn.Omega(a)**2

def Gamma_tot(a,m_core,m_atm,t):
    return -(1.36+.62*ct.beta_sigma+.43*ct.beta_t)*Gamma_0(a,m_core,m_atm,t)

#inputs in AU, Mearth, Myr. Note that the derivatives should be in these units.
def a_dot(a,m_core,m_core_dot,m_atm,m_atm_dot,t): #a note: i could combine the m_core/m_atm and m_core_dot/m_atm_dot parameters
    
    #return 0 #uncomment for a no-migration run.
    
    m_core=m_core*ct.Mearth
    m_atm=m_atm*ct.Mearth

    if(check_gap(a,(m_core+m_atm))):
        return 0 #placeholder for now. need to implement viscous migration rate.

    cgs_a_dot= 2*(Gamma_tot(a,m_core,m_atm,t)-(m_atm_dot+m_core_dot)*ct.Mearth/ct.Myr_to_s*(ct.G*ct.Msun)**.5*(a*ct.AU_to_cm)**.5)/((m_core+m_atm)*(ct.G*ct.Msun)**.5*(a*ct.AU_to_cm)**-.5)
    return cgs_a_dot*ct.Myr_to_s/ct.AU_to_cm
    
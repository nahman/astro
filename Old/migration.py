import constants as ct
import basic_equations as eqn 

'''
to do: 
-make discontinous transitions of a_dot smoother.
'''

#viscous migration speed of planet. (Type II). in AU/Myr
a_dot_visc=-ct.alpha_ss*(ct.k/(ct.mu*ct.m_H)*260*(ct.AU_to_cm)**.5/(ct.G*ct.Msun)**.5)*ct.Myr_to_s/ct.AU_to_cm

#normalization torque (inputs in AU, GRAMS, myr)
def Gamma_0(a,m_core,m_atm,t): 
    m = m_core+m_atm
    return (m/ct.Msun)**2*(eqn.H(a)/(a*ct.AU_to_cm))**-2*eqn.Sigma_gas_t(a,t)*(a*ct.AU_to_cm)**4*eqn.Omega(a)**2

#total torque
def Gamma_tot(a,m_core,m_atm,t):
    return -(1.36+.62*ct.beta_sigma+.43*ct.beta_t)*Gamma_0(a,m_core,m_atm,t)

#time derivative of planet's orbital distance. inputs in AU, Mearth, Myr. 
def a_dot(a,m_core,m_atm,t): 
    if ct.const_mig:
        return 0 

    m_core=m_core*ct.Mearth
    m_atm=m_atm*ct.Mearth
    
    if(eqn.check_gap(a,(m_core+m_atm))):
        return 0 #turning of viscous migration
        return a_dot_visc

    cgs_a_dot= 2*Gamma_tot(a,m_core,m_atm,t)/((m_core+m_atm)*(ct.G*ct.Msun)**.5*(a*ct.AU_to_cm)**-.5)
    return cgs_a_dot*ct.Myr_to_s/ct.AU_to_cm
    

#Eve's a_dot expression
def a_dot_alt(a,m_core,m_atm,t):
    if(eqn.check_gap(a,(m_core+m_atm)*ct.Mearth)):
        print('gap formed')
        return 0
    dadt = -(1.36+0.62*ct.beta_sigma+0.43*ct.beta_t)*2*ct.G**2*((m_atm+m_core)*ct.Mearth*eqn.Sigma_gas_t(a,t)/eqn.c_s(a)**3)*(eqn.H(a)/(a*ct.AU_to_cm))
    dadt = dadt*ct.Myr_to_s/ct.AU_to_cm
    return dadt


'''
testing
'''

#print(Gamma_0(1,2,1,1e-1))
print(a_dot(1,2,1,1e-1))


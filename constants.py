import numpy as np


'''
Constants that should never change
'''
#Physical
G=6.67e-8
Msun=2e+33
Mearth=6e+27
k=1.38e-16
m_H=1.67e-24
rho_solid=5.0
mu=2.374 #mean molecular weight in hydrogen masses
mu_rcb=2.0 #I'm not sure if this should be 2 or 2.374
T_rcb = 2500.0 
grad_adb=.17 #adiabatic temperature gradient at rcb
m_jup=300 #in mearths

#Conversions
AU_to_cm=1.5e+13
Myr_to_s=3.1557e+13

#Math
pi=np.pi

'''
Parameters that might be changed
'''

tau_fr=1.0 #dimensionless friction time 
C_acc = 1.0
gap_param=1.0 #If sigma_gap/sigma_0 is less than this, we say a gap exists, and core accretion stops.
fudge=2.0 #fudge factor for GCR equations. FC14 states f is around 2 for 1AU to 5AU, 3 for less than 1AU
Z=2.0e-2 #gas metallicity 

beta_sigma=1.5 #these are parameters for calculating the total torque
beta_t=.429

Pcol_fudge=3 #divisor to Pcol

core_smooth=10 #higher numbers=rounder points
atm_smooth=70

'''
Gas parameters
'''

depletion_time=5 #in Myrs. how long does the gas disk last.
depletion_factor=1
alpha_ss=1.0e-4 #turbulent parameter

'''
Run parameters
'''
Sigma_pow=True #enable for time dependent power scaling of sigma_gas. otherwise, efolding.

const_core=False
const_atm=False
const_mig=False
const_gas=False

#start/stop times, in log space
t_init = -3
t_end= 1

step=10000

#initial conds, in Mearth, AU
y0=[.01,0,10] 


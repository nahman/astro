import numpy as np

'''
Constants that should never change
'''
#Physical Constants, cgs
G=6.67e-8 #gravitational constant
Msun=2e+33
Mearth=6e+27
k=1.38e-16 #boltzmann constant
sb=5.67e-5 #stefan boltzmann constant
m_H=1.67e-24 #hydrogen mass
mu=2.37 #mean molecular weight in hydrogen masses
mu_rcb=2.37 #I'm not sure if this should be 2 or 2.374
T_rcb = 2500.0 #temperature of atm at rcb
grad_adb=.17 #adiabatic temperature gradient at rcb
grad_adb_clean=.25
m_jup=300 #Jupiter mass in mearths

#Conversions
AU_to_cm=1.5e+13 #one AU is _ cm
Myr_to_s=3.1557e+13 #one Myr is _ s

#Math
pi=np.pi

#Migrational Torque
beta_sigma=1.0
beta_t=.429



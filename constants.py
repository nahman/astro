import numpy as np
import itertools as it

iterator = it.count()
index=iterator.__next__()

def get_next():
    global index
    index=iterator.__next__()
    return


'''
Constants that should never change
'''
#Physical Constants, cgs
G=6.67e-8 #gravitational constant
Msun=2e+33
Mearth=6e+27
k=1.38e-16 #boltzmann constant
m_H=1.67e-24 #hydrogen mass
mu=2.374 #mean molecular weight in hydrogen masses
mu_rcb=2.0 #I'm not sure if this should be 2 or 2.374
T_rcb = 2500.0 #temperature of atm at rcb
grad_adb=.17 #adiabatic temperature gradient at rcb
m_jup=300 #Jupiter mass in mearths

#Conversions
AU_to_cm=1.5e+13 #one AU is _ cm
Myr_to_s=3.1557e+13 #one Myr is _ s

#Math
pi=np.pi

#Gap Opening
gap_param = 1.0 #NOT IN USE RN If sigma_gap/sigma_0 is less than this, we say a gap exists, and core accretion stops.

#Migrational Torque
beta_sigma=1.5 #these are parameters for the torque calculation in migration.py
beta_t=.429


'''
todo: 

-implement a mixed power - exp gas depletion. Once Sigma drops by 2 OOM, power becomes exp.
-test params (relevant values listed in parentheses above)
-Case 1: Once M_core = M_Iso, switch to type II migration. Vary alpha logspace (-2,-3,-4,-5) 
-make a m_planet, period plot. Then overlay the gap opening criterion line, make sure it intersects before type II migration starts.
-Case 2: Try this case: Always use type I until M_p \sim M_jup
-What parameters cause core growth then migration, vs the opposite?
-Case 3: Keep Type I till Jupiter formation, then use type II.


In all cases we want a m_p, period plot.
For each case we should try both scalings.

'''

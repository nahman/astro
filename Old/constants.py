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
rho_solid=5.0 #density of core stuff
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

'''
Parameters that might be changed. []
'''

tau_fr = 1.0 #dimensionless friction time 
C_acc = 1.0 # mass accretion rate parameter in OK12: set to 1 for single-body growth.
gap_param = 1.0 #If sigma_gap/sigma_0 is less than this, we say a gap exists, and core accretion stops.
fudge = 2.0 #fudge factor for GCR equations. FC14 states f is around 2 for 1AU to 5AU, 3 for less than 1AU
Z = 2.0e-2 #gas metallicity 

beta_sigma=1.5 #these are parameters for the torque calculation in migration.py
beta_t=.429

Pcol_fudge=3 #divisor to Pcol. Setting to 3 matches this Pcol with Eve's (approximately)


core_smooth=1e+10 #These params control how quickly the time derivatives of core mass and atm mass decay to 0 after the core/atm stops growing.
atm_smooth=1e+10  #The decay is exponential in form, and mainly was implemented to keep the derivatives continous.
                  #higher numbers=sharper points. be careful with low numbers; they tend to inflate the final value of the parameter.

iso_mass=20 #core isolation mass in mearths

'''
Gas parameters
'''

depletion_time=5   #in Myrs. How long does the gas disk last. (.1,1,10)
depletion_factor=.1 #multiplied against the gas surf density to raise/lower it. 
alpha_ss=1e-4   #turbulent parameter


'''
Run parameters
'''
reg='auto' #manually set the regime we are in when we calculate Pcol
          #values: set, hyp, reg. set to auto for automatic selection (via the max function).

Sigma_pow=True #enable for time dependent power scaling of sigma_gas. otherwise, the code uses an efolding timescale.

#enable for the respective parameter to be held constant during the run.
const_core=False 
const_atm=False
const_mig=False
const_gas=False

#start/stop times, in log space, ie 10^t_init to 10^t_end Myrs.
t_init = -3  #params to test (-3,-2,-1)
t_end= 1

#number of steps for odeint to use
step=1000

#initial conds [coremass, atmass, orbital dist] in Mearth, Mearth, AU
y0=[1e-2,0,1] 

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

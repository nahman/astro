import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

import atmass as am
import basic_equations as eqn
import constants as ct
import coremass as cm
import migration as mig


new_sol=[]


def diff_system(y,t):

    m_core=y[0]
    m_atm=y[1]
    a=max(1e-2,y[2]) #prevent a from becoming negative
    print('a ' + str(a))
    print('m_core ' + str(m_core))
    print('m_atm '+str(m_atm))

    
    m_core_dot = cm.M_core_dot(a,m_core,m_atm,t) 
    m_atm_dot = am.m_atm_dot_dusty(a,m_core,m_core_dot,m_atm,t)
    a_dot = mig.a_dot(a,m_core,m_atm,t)
    
    dydt = [m_core_dot,m_atm_dot,a_dot]
    return dydt


def graph(sol,t,names,log): #I only have 6 colors, so max 6 thingies
    colors = ['b','g','r','c','m','y']

    for i in range(len(names)):
        plt.plot(t, sol[:,i], colors[i], label=names[i])

    plt.legend(loc='best')
    plt.grid()
    if log:
        plt.xscale('log')
        plt.yscale('log')
    plt.show()


#initial conds
y0=[.01,0,10]

#solving times
t_log=np.logspace(-3,1,100000) 

sol = odeint(diff_system,y0,t_log)

graph(sol,t_log,['m_core','m_atm','a'],True)

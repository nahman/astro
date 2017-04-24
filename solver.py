import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

import atmass as am
import basic_equations as eqn
import constants as ct
import coremass as cm
import migration as mig

def diff_system(y,t):

    m_core=y[0]
    m_atm=y[1]
    a=max(1e-3,y[2]) #prevent a from becoming negative

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
    
    plt.axvline(x= cm.t_gap) #marks when gap forms

    plt.show()

def solve_graph(t_init, t_end, y, y0, step, names,log=True): #linear not implemented yet
    if log:
        t = np.logspace(t_init,t_end,step)
    sol = odeint(y,y0,t)
    graph(sol,t,names,log)

'''
for x in range(len(sol)):
    print(sol[x],t_log[x])
'''
#print('final values')
#print(sol[-1])

solve_graph(ct.t_init,ct.t_end,diff_system,ct.y0,ct.step,['m_core','m_atm','a'])


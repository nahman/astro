import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

import atmass as am
import basic_equations as eqn
import constants as ct
import coremass as cm
import migration as mig

#The system of diffeqs we are solving. Takes m_core, m_atm, a, and returns their derivatives.
def diff_system(y,t):

    m_core=y[0]
    m_atm=y[1]
    a=max(1e-3,y[2]) #prevent a from becoming negative

    m_core_dot = cm.M_core_dot(a,m_core,m_atm,t) 
    m_atm_dot = am.m_atm_dot_dusty(a,m_core,m_core_dot,m_atm,t)
    a_dot = mig.a_dot(a,m_core,m_atm,t)
    
    dydt = [m_core_dot,m_atm_dot,a_dot]
    return dydt

#Uses matplotlib to graph any generic solution outputted by odeint.
def graph(sol,t,names,log): 
    colors = ['b','g','r','c','m','y']

    for i in range(len(names)):
        plt.plot(t, sol[:,i], colors[i%6], label=names[i])

    plt.legend(loc='best')
    plt.grid()
    if log:
        plt.xscale('log')
        plt.yscale('log')
    plt.xlabel('time, Myr')
    plt.ylabel('core mass, Mearth')
    
    #plt.axvline(x= cm.t_gap) #marks when gap forms

    plt.show()

#Function that takes care of both solving diff_system with odeint and graphing to solution with matplotlib
#takes, start/stop times, the diff system in question, initial conditions, the desired of number of steps, 
#and the names of each equation you are trying to solve for. The last paramter decides whether you plot in
#log space or linear space.

def solve_graph(t_init, t_end, y, y0, step, names,log=True): #linear not implemented yet
    if log:
        t = np.logspace(t_init,t_end,step)
    sol = odeint(y,y0,t)
    graph(sol,t,names,log)

    #get the value of the solution to the diff eq system at t_end
    print('final values')
    print(sol[-1])



solve_graph(ct.t_init,ct.t_end,diff_system,ct.y0,ct.step,['m_core','m_atm','a']) #I should wrap this in some main function?


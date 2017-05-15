import numpy as np
import constants as ct
import disk as d 
import planet as pl
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def diff_wrapper(planet:pl.Planet):
    
    def diff_system(y,t):
        nonlocal planet
        planet.redefine(y[0],y[1],y[2])
        m_core_dot = planet.m_core_dot(t)
        m_atm_dot = planet.m_atm_dot_dusty(t)
        a_dot = planet.a_dot(t)
        print(planet.m_atm)
        if planet.a<=1e-2:
            planet.a=1e-2
            a_dot=0

        return [m_core_dot,m_atm_dot,a_dot]

    return diff_system
    

def graph(t,sol,names,index=-1,name=-1,color=-1): 
    colors = ['b','g','r','c','m','y']
    
    if index==-1:
        for i in range(len(names)):
            plt.plot(t, sol[:,i], colors[i%6], label=names[i])
        return 
    plt.plot(t,sol[:,index],colors[color],label=name)

def format_graph(title=None,xlabel=None,ylabel=None,showline=False,log=True):
    plt.legend(loc='best')
    plt.grid()
    if log:
        plt.xscale('log')
        plt.yscale('log')
    if xlabel !=None:
        plt.xlabel(xlabel)
    if ylabel !=None:
        plt.ylabel(ylabel)
    if title !=None:
        plt.title(title)
    if showline:
        plt.axvline(x= cm.t_stop) #marks when gap forms
    plt.tight_layout()
    plt.show()



y0=[1e-2,0,1]
t=np.logspace(-3,3)

disk=d.Disk('pow',1,1e-4,t_init=1e-3)
planet =pl.Planet(y0[0],y0[1],y0[2],10,disk)


sol = odeint(diff_wrapper(planet),planet.state(),t)

graph(t,sol,['m_c','m_a','a'])
format_graph()
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

    if eqn.check_gap(a,(m_atm+m_core)*ct.Mearth):
        m_atm=ct.m_jup-m_core

    #print(m_core,m_atm,a)
    m_core_dot = cm.M_core_dot_iso(a,m_core,m_atm,t) 
    m_atm_dot = am.m_atm_dot_dusty(a,m_core,m_core_dot,m_atm,t)
    a_dot = mig.a_dot(a,m_core,m_atm,t)
    dydt = [m_core_dot,m_atm_dot,a_dot]
    return dydt

#Uses matplotlib to graph any generic solution outputted by odeint.
#pass in index to select a particular part of sol to graph.

def graph(t,sol,names,index=-1,name=-1,color=-1): 
    colors = ['b','g','r','c','m','y']
    
    if index==-1:
        for i in range(len(names)):
            plt.plot(t, sol[:,i], colors[i%6], label=names[i])
        return 
    plt.plot(t,sol[:,index],colors[color],label=name)
    
    #plt.plot(t,np.vectorize(eqn.Sigma_factor)(t),'m',label='Sigma_gas factor')
    
#Run this at the end to show the graph
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


#Function that takes care of both solving diff_system with odeint and graphing to solution with matplotlib
#takes, start/stop times, the diff system in question, initial conditions, the desired of number of steps, 
#and the names of each equation you are trying to solve for. The last paramter decides whether you plot in
#log space or linear space.

def solve(t_init, t_end, y, y0, step,log=True): #linear not implemented yet

    if log:
        t = np.logspace(t_init,t_end,step)
    else:
        raise Exception('I forgot to implement the linear case')
    
    sol = odeint(y,y0,t)
    #get the value of the solution to the diff eq system at t_end
    print('final values')
    print(sol[-1])
    return (t, sol)

#Generator that varies a single parameter and yields the solutions one at a time

def loop(vals,param_title):
    i=0
    for x in vals:
        ct.iso_mass=x  #change this to whatever
        sol = solve(ct.t_init,ct.t_end,diff_system,ct.y0,ct.step)
        i+=1
        yield sol
    return None


#graphs one of the solution functions (core, atm, a) over a number of parameter values. 'single panel'

def sp_graph(vals,param_title,index):
    sol_gen=loop(vals,param_title)
    
    sol=next(sol_gen)
    i=0
    while True:
        plt.loglog(sol[0],sol[1][:,index],label=param_title+' '+str(vals[i]))
        i+=1
        try:
            sol=next(sol_gen)
        except Exception as e:
            print(e)
            break
    


#Similar to sp_graph, but plots all three solution functions in a three panel plot.
def mp_graph(vals, param_title):
    
    sol_dict={}
    
    sol_gen = loop(vals,'Iso Mass')
    sol=next(sol_gen)

    fig = plt.figure()

    i=0
    while True:
        plt.subplot(311)
        name = 'm_core '+param_title+' '+str(vals[i])
        sol_dict[name], = plt.loglog(sol[0],sol[1][:,0],label='IsoM '+str(vals[i]))
        
        plt.subplot(312)
        name = 'm_atm '+param_title+' '+str(vals[i])
        sol_dict[name], = plt.loglog(sol[0],sol[1][:,1],label='IsoM '+str(vals[i]))
        
        plt.subplot(313)
        name = 'a '+param_title+' '+str(vals[i])
        sol_dict[name], = plt.loglog(sol[0],sol[1][:,2],label='IsoM '+str(vals[i]))
        
        i+=1
        
        try:
            sol=next(sol_gen)
        except Exception as e:
            print(e)
            break
        
    
    plt.subplot(311)
    format_graph(ylabel='Mass (Mearth)',title='M_core')
    plt.legend(handles=[sol_dict['m_core '+param_title+' '+str(x)] for x in vals],loc='upper right')
    
    plt.subplot(312)
    format_graph(ylabel='Mass (Mearth)',title='M_atm')
    plt.legend(handles=[sol_dict['m_atm '+param_title+' '+str(x)] for x in vals],loc='right')
    
    plt.subplot(313)
    format_graph(ylabel='Dist (AU)',title='Orbital Radius')
    plt.legend(handles=[sol_dict['a '+param_title+' '+str(x)] for x in vals],loc='right')

    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    plt.xlabel('Time, Myrs')

'''
vals=[5,10,15,20,25]

mp_graph(vals,'Iso M')
plt.suptitle('Effect of Iso M on Planet Evolution')
plt.show()
'''
t = np.logspace(ct.t_init,ct.t_end,ct.step)
sol = solve(ct.t_init,ct.t_end,diff_system,ct.y0,ct.step)
print(sol)
graph(sol[0],sol[1],['mc','ma','a'])

format_graph()
plt.show()
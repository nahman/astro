import numpy as np
import constants as ct 
import basic_equations as eqn 
import atmass as am
import coremass as cm
import migration as mig
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def diff_system(y,t):
    m_core, m_atm, a = y
    #print(m_core,m_atm,a,t)
    m_core_dot = cm.M_core_dot(a,m_core,m_atm)
    m_atm_dot = am.m_atm_dot_dusty(a,m_core,m_core_dot,m_atm,t)
    a_dot = mig.a_dot(a,m_core,m_core_dot,m_atm,m_atm_dot,t)
    dydt = [m_core_dot,m_atm_dot,a_dot]
    return dydt

y0=[1.0,.001,1.0]

t0=np.linspace(.00001,.0002,1000)

sol = odeint(diff_system,y0,t0)
print(sol)

plt.plot(t0, sol[:,1], 'b', label='m_atm')
plt.plot(t0, sol[:,0], 'g', label='m_core')
plt.plot(t0, sol[:,2], 'r', label='a')
plt.legend(loc='best')
plt.grid()
#plt.xscale('log')
plt.show()
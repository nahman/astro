import numpy as np 
import matplotlib.pyplot as plt
import constants as ct


AU_to_cm=1.5e+13 #one AU is _ cm
Myr_to_s=3.1557e+13 #one Myr is _ s

R1=1
C=52000
gamma=1
alpha = 1e-2

def temp(a):                                #temp of disk. takes AU, K
    return 260*(a**-.5)

def c_s(a):                                 #sound speed. input AU, output CGS
    return (ct.k*temp(a)/(ct.mu*ct.m_H))**.5

def Omega(a):                               #Keplerian angular velocity. input in AU, output rad/s
    return (ct.G*ct.Msun/(a*ct.AU_to_cm)**3)**.5

def H(a):                                   #disk scale height, input in AU, output in cgs
    return c_s(a)/Omega(a)

def nu_alt(r): #input AU, output in AU^2/Myr
    return alpha*c_s(r)*H(r)*AU_to_cm**-2*Myr_to_s

'''
def nu(r):
    return r**gamma
'''

nu1=nu_alt(R1)


def r(R):
    return R/R1

t_s = 1/(3*(2-gamma)**2)*R1**2/nu1

def T(t):
    return t/t_s+1


def Sigma(R,t):

    ret = C/(3*np.pi*nu1*r(R)**gamma) * T(t)**(-(5/2-gamma)/(2-gamma))*np.exp(-r(R)**(2-gamma)/T(t))
    return ret

print(10**2/nu_alt(10))


R=np.logspace(0,4)

t_vals=np.arange(0,10,1)

for t in t_vals:
    plt.loglog(R,np.vectorize(lambda R: Sigma(R,t))(R),label='t='+str(t)+'Myr')
plt.xlabel('distance, AU')
plt.ylabel('surface density, g/cm^2')
plt.title('Sigma_gas at various times')
#plt.loglog(R,100*(R/10)**-1.5,label='-1.5 power law')


plt.legend(loc='best')
plt.show()



'''

times = np.logspace(-3,1)
plt.loglog(times,np.vectorize(lambda t: Sigma(1,t))(times),label='new sigma_gas')
plt.loglog(times,np.vectorize(lambda t: .55*t**(-5/4))(times),label='power depletion')
'''

'''
R=np.logspace(0,2)
R1_vals = [1,2,5,10,100]
t_val=10

c_list=[]

for val in R1_vals:
    R1=val
    C=1
    c_list.append(3173/Sigma(1,0))


for i in range(len(R1_vals)):
    R1=R1_vals[i]
    C=c_list[i]

    plt.loglog(R,np.vectorize(Sigma)(R,t_val),label='R1='+str(R1)+'AU')

plt.xlabel('distance, AU')
plt.ylabel('surface density, g/cm^2')
plt.title('Sigma_gas at t='+str(t_val)+'myr over various R1s')

plt.ylim(1e-10,1e5)

plt.legend(loc='best')
plt.show()

'''
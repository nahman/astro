import numpy as np
import constants as ct
import basic_equations as eqn


'''
Expressions for Pcol in different regimes
'''

tau_fr=ct.tau_fr

#Condition checking for what regime we are in:

def check_reg(a,m):
    zeta=eqn.zeta_w(a,m)
    reg=''
    if tau_fr<=min(12/zeta**3,2):
        reg='set'
    elif tau_fr>=zeta:
        reg='3b'
    else:
        reg='hyp'
    return (reg,zeta) #returning zeta so I don't have to call the function again. Probably a more elegant way to do this.


#Condition checking for gap opening (from Eugene's paper)

def check_gap(a,m):
    q = m/ct.Msun

    h_r= eqn.H(a)/(a*ct.AU_to_cm) #disk aspect ratio

    if 1e-4<=q<5e-3:

        if .14*(q/1e-3)**-2.16*(ct.alpha_ss/1e-2)*(h_r/.05)**6.61 <= ct.gap_param:
            return True
        return False

    elif 5e-3<=q<=1e-2:

        if (q/5e-3)**-1.00*(ct.alpha_ss/1e-2)**1.26*(h_r/.05)**5.12 <= ct.gap_param:
            return True
        return False
    elif q<1e-4:
        return False #I assume low mass planets don't form gaps
    else:

        raise Exception('check_gap failed. the mass of the planet exceeds upper bounds set for the calculation of gap formation')



#Pcol function. This one actually checks regimes.

def Pcol_alt(a,m): #this stuff is all in normalized hill units #inputs in AU, grams
    reg,zeta = check_reg(a,m)
    alpha=eqn.alpha_core(a,m)
    b=None   #impact param
    v_a=None #approach velocity
    if reg=='set':
        coeff = [1,2/3*zeta,0,-8*tau_fr]
        roots = np.roots(coeff)
        b=roots.real[abs(roots.imag)<1e-5][0]
        v_a=3/2*b+zeta
        '''
        exponential smoothing factor for b, unused for now
        '''
        f = np.exp(-(tau_fr/min(12/zeta**3,2))**.65)
    elif reg=='hyp':
        v_a=zeta*(1+4*tau_fr**2)**.5/(1+tau_fr)**2
        b=alpha*(1+6/(alpha*v_a**2))**.5
        '''
        no smoothing factor provided OK12
        '''
    else:
        b=1.7*alpha+1/tau_fr
        v_a=3.2
        '''
        exponential smoothing factor for b, unused for now
        '''
        f=np.exp(-(.7*zeta/tau_fr)**5)
    return 2*b*eqn.r_hill(a,m)*v_a*eqn.v_hill(a,m)

#Alternative Pcol function. This one takes the max (which frankly seems to be wasteful in terms of computation.)

def Pcol(a,m): #input in AU, grams

    zeta = eqn.zeta_w(a,m)
    r_H = eqn.r_hill(a,m)
    v_H = eqn.v_hill(a,m)
    alpha=eqn.alpha_core(a,m)

    #settling calc
    coeff = [1,2/3*zeta,0,-8*tau_fr]
    coeff=np.nan_to_num(coeff)
   # print(coeff)
    roots = np.roots(coeff)
    b_set=roots.real[abs(roots.imag)<1e-5][0]
    v_a_set=3/2*b_set+zeta
    f_set = np.exp(-(tau_fr/min(12/zeta**3,2))**.65) #smoothing factor

    Pcol_set = 2*b_set*f_set*v_a_set 
    #print(Pcol_set)

    #hyperbolic calc
    v_a_hyp=zeta*(1+4*tau_fr**2)**.5/(1+tau_fr)**2
    b_hyp=alpha*(1+6/(alpha*v_a_hyp**2))**.5

    P_col_hyp = 2*b_hyp*v_a_hyp

    #three body calc
    b_3b=1.7*alpha+1/tau_fr
    v_a_3b=3.2
    f_3b=np.exp(-(.7*zeta/tau_fr)**5) #smoothing factor
    
    P_col_3b = 2*b_3b*f_3b*v_a_3b

    return max(Pcol_set,P_col_hyp,P_col_3b)
    

#Expression for derivative

def M_core_dot(a,m_core,m_atm): #inputs in M_earth, AU
    m=(m_core+m_atm)*ct.Mearth #planet mass in cgs

    if check_gap(a,m)==False:

        return ct.C_acc*Pcol(a,m)*eqn.Sigma_mmsn(a)*eqn.r_hill(a,m)*eqn.v_hill(a,m)*ct.Myr_to_s/ct.Mearth

    return 0 #gap formed, no more accretion

'''
testing
'''
#print (M_core_dot(1,1,0))
#print (M_core_dot(1,1,1))
#print (M_core_dot(1,100,0))


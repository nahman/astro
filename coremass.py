import numpy as np
import constants as ct
import basic_equations as eqn

tau_fr=ct.tau_fr
t_gap=0 #gets set to the time at which the planet opens a gap. Used later during plotting.



#These functions are all from OK12, Appendix A

'''
Regime Checking
'''

#What regime is the planet forming in? Settling, Hyperbolic, or Three-body? 
#Input in AU, grams. Returns zeta as well as the regime to save a calculation.
def check_reg(a,m):
    zeta=eqn.zeta_w(a,m)
    reg=''
    if tau_fr<=min(12/zeta**3,2):
        reg='set'
    elif tau_fr>=zeta:
        reg='3b'
    else:
        reg='hyp'
    return (reg,zeta) 


'''
Expressions for Pcol
'''

#Pcol function. This one actually checks regimes.

def Pcol_alt(a,m):                      #inputs in AU, grams. outputs in hill units.
    reg,zeta = check_reg(a,m)
    alpha=eqn.alpha_core(a,m)
    b=None                              #impact param
    v_a=None                            #approach velocity
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


#Alternative Pcol function (which is being used currently). This one takes the max 
#(which frankly seems to be wasteful in terms of computation).

def Pcol(a,m): #input in AU, grams. outputs in hill units.

    zeta = eqn.zeta_w(a,m)
    r_H = eqn.r_hill(a,m)
    v_H = eqn.v_hill(a,m)
    alpha=eqn.alpha_core(a,m)

    #settling calc
    coeff = [1,2/3*zeta,0,-8*tau_fr]
    coeff=np.nan_to_num(coeff)
    roots = np.roots(coeff)
    b_set=roots.real[abs(roots.imag)<1e-5][0]
    v_a_set=3/2*b_set+zeta
    f_set = np.exp(-(tau_fr/min(12/zeta**3,2))**.65) #smoothing factor

    Pcol_set = 2*b_set*f_set*v_a_set 

    #hyperbolic calc
    v_a_hyp=zeta*(1+4*tau_fr**2)**.5/(1+tau_fr)**2
    b_hyp=alpha*(1+6/(alpha*v_a_hyp**2))**.5

    P_col_hyp = 2*b_hyp*v_a_hyp

    #three body calc
    b_3b=1.7*alpha+1/tau_fr
    v_a_3b=3.2
    f_3b=np.exp(-(.7*zeta/tau_fr)**5) #smoothing factor
    
    P_col_3b = 2*b_3b*f_3b*v_a_3b

    #check if we should force a regime or let the max do its thing.
    if ct.reg=='auto':
        return max(Pcol_set,P_col_hyp,P_col_3b)
    elif ct.reg=='set':
        return Pcol_set
    elif ct.reg=='hyp':
        return P_col_hyp
    elif ct.reg=='3b':
        return P_col_3b
    else:
        raise Exception('error in Pcol. no regime with such name')
    

#Expression for derivative. From OK12, eqn 1.

def M_core_dot(a,m_core,m_atm,t): #inputs in M_earth, AU
    
    if ct.const_core:
        return 0
   
    m=(m_core+m_atm)*ct.Mearth #planet mass in cgs

    #expression for m_core_dot. in mearth/myr.
    m_dot = ct.C_acc*Pcol(a,m)/ct.Pcol_fudge*eqn.Sigma_mmsn(a)*eqn.r_hill(a,m)*eqn.v_hill(a,m)*ct.Myr_to_s/ct.Mearth 

    #Check that the planet has not formed a gap. Gap formation halts accretion.
    if eqn.check_gap(a,m)==False:
        
        return m_dot
    
    #Otherwise, record the time of gap_opening.
    global t_gap
    
    if t_gap==0:
        t_gap=t
    
    #gap formed, no more accretion. Thus m_dot decays to 0.
    #This decay is caused by the exponential factor.
    return m_dot*np.exp(-ct.core_smooth*(t-t_gap)) 



'''
testing
'''
#print (M_core_dot(1,1,0))
#print (M_core_dot(1,1,1))
#print (M_core_dot(1,100,0))
#print(Pcol(1,ct.Mearth))


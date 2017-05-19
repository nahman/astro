#A class that makes the identification of analyses easier. 
import json

def sci(num):
    return "{:.1e}".format(num)

'''
A Potential Problem
Planet_ID's dont depend on initial conds. This is an easy (but tedious fix). Not sure if it's worth the effort as of now.
'''

class Planet_ID: 

    def __init__(self,tau_fr=None,alpha=None,m_iso=None,mode=None,t_init=None,_dict=None):
        #more parameters can be added here. wonder if there is a more concise/pythonic way...
        if _dict!=None:
            alpha= _dict['alpha']
            m_iso= _dict['m_iso']
            mode= _dict['mode']
            tau_fr= _dict['tau_fr']
            t_init= _dict['t_init']

        #below is the standard non-dict initialization
        self.alpha=alpha
        self.m_iso=m_iso
        self.mode=mode
        self.tau_fr=tau_fr
        self.t_init=t_init
    

    def name(self):
        return 'alpha: '+str(sci(self.alpha))+' | m_iso: '+str(self.m_iso)+' | mode: '+str(self.mode)+' | tau: '+str(self.tau_fr)+' | t_init: '+str(sci(self.t_init))

    def __str__(self):
        params = {}
        params['tau_fr']=self.tau_fr    
        params['alpha']=self.alpha
        params['m_iso']=self.m_iso
        params['mode']=self.mode
        params['t_init']=self.t_init
        return json.dumps(params)
    
    def fromstr(self,_str):
        d = json.loads(_str)
        self.tau_fr=d['tau_fr']
        self.alpha=d['alpha']
        self.m_iso=d['m_iso']
        self.mode=d['mode']
        self.t_init=d['t_init']

    def __hash__(self):
        return hash((self.alpha,self.m_iso,self.mode,self.tau_fr,self.t_init))

    def __eq__(self,other):
        return (self.alpha,self.m_iso,self.mode,self.tau_fr,self.t_init)==(other.alpha,other.m_iso,other.mode,other.tau_fr,other.t_init)


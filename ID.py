#A class that makes the identification of analyses easier. 
import json

def sci(num):
    return "{:.1e}".format(num)

class Planet_ID: 

    def __init__(self,tau_fr=None,alpha=None,m_iso_mode=None,mode=None,t_init=None,t_dep=None,atmos:'clean or dusty'=None,case=None,_dict=None):
        #more parameters can be added here. wonder if there is a more concise/pythonic way...
        if _dict!=None:
            alpha= _dict['alpha']
            m_iso_mode= _dict['m_iso_mode']
            mode= _dict['mode']
            tau_fr= _dict['tau_fr']
            t_init= _dict['t_init']
            t_dep=_dict['t_dep']
            atmos=_dict['atmos']
            case=_dict['case']

        #below is the standard non-dict initialization
        self.alpha=alpha
        self.m_iso_mode=m_iso_mode
        self.mode=mode
        self.tau_fr=tau_fr
        self.t_init=t_init
        self.t_dep=t_dep
        self.atmos = atmos
        self.case=case
    

    def param_list(self):
        t_dep_or_init='\nt_dep: '+str(self.t_dep)
        if self.mode!='exp':
            t_dep_or_init='\nt_init: '+str(sci(self.t_init))

        ret =('alpha: '+ str(sci(self.alpha)) +
        '\nm_iso_mode: '+str(self.m_iso_mode) +
        '\ndisk_mode: '+str(self.mode) +
        '\ntau: '+str(self.tau_fr) +
        t_dep_or_init +
        '\natmos: '+str(self.atmos) +
        '\ncase: '+str(self.case))
           
        return ret

    def file_name(self):
        t_dep_or_init=' t_dep '+str(self.t_dep)
        if self.mode!='exp':
            t_dep_or_init=' t_init '+str(sci(self.t_init))

        return ('alpha '+str(sci(self.alpha))+
        ' m_iso_mode '+str(self.m_iso_mode)+
        ' disk_mode '+str(self.mode)+
        ' tau '+str(self.tau_fr)+
        t_dep_or_init+
        ' atmos '+str(self.atmos)+
        ' case '+str(self.case))
        

    def __str__(self):
        params = {}
        params['tau_fr']=self.tau_fr    
        params['alpha']=self.alpha
        params['m_iso_mode']=self.m_iso_mode
        params['mode']=self.mode
        params['t_init']=self.t_init
        params['t_dep']=self.t_dep
        params['atmos']=self.atmos
        params['case']=self.case
        return json.dumps(params)
    
    def fromstr(self,_str):
        d = json.loads(_str)
        self.tau_fr=d['tau_fr']
        self.alpha=d['alpha']
        self.m_iso_mode=d['m_iso_mode']
        self.mode=d['mode']
        self.t_init=d['t_init']
        self.t_dep=d['t_dep']
        self.atmos=d['atmos']
        self.case=d['case']

    def __hash__(self):
        return hash((self.alpha,self.m_iso_mode,self.mode,self.tau_fr,self.t_init,self.t_dep,self.atmos,self.case))

    def __eq__(self,other):
        return (self.alpha,self.m_iso_mode,self.mode,self.tau_fr,self.t_init,self.atmos,self.case) \
        == (other.alpha,other.m_iso,other.mode,other.tau_fr,other.t_init,other.t_dep,other.atmos,other.case)


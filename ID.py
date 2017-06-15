#A class that makes the identification of analyses easier. 
import json

def sci(num):
    num = float(num)
    return "{:.1e}".format(num)

class Planet_ID: 

    def __init__(self,
                tau_fr : float = None,
                alpha : float = None,
                m_iso_mode=None,
                t_s=None,
                atmos:'clean or dusty'=None,
                case=None,
                _dict:'pass in a dictionary containing all relevant parameter values for instant init'=None
                ):

        #used with optional dictionary normalization
        if _dict!=None:
            alpha= _dict['alpha']
            m_iso_mode= _dict['m_iso_mode']
            tau_fr= _dict['tau_fr']
            t_s= _dict['t_s']
            atmos=_dict['atmos']
            case=_dict['case']

        #below is the standard non-dict initialization
        self.alpha=sci(float(alpha))
        self.m_iso_mode=m_iso_mode
        self.tau_fr=sci(float(tau_fr))
        self.t_s=t_s
        self.atmos = atmos
        self.case=case
    
    #outputs the parameter values into a broken string that can be added to plots
    def param_list(self):

        ret =('alpha: '+ str(sci(self.alpha)) +
        '\ntau: '+str(sci(self.tau_fr)) +
        '\nt_s: '+str(self.t_s) +
        '\nm_iso_mode: '+str(self.m_iso_mode) +
        '\natmos: '+str(self.atmos) +
        '\ncase: '+str(self.case))
           
        return ret

    #outputs the parameter values into a string that can be used for naming
    def file_name(self):

        return ('alpha '+str(sci(self.alpha))+
        ' tau '+str(sci(self.tau_fr))+
        ' t_s '+str(self.t_s)+
        ' m_iso_mode '+str(self.m_iso_mode)+
        ' atmos '+str(self.atmos)+
        ' case '+str(self.case))
        
    def __str__(self):
        params = {}
        params['tau_fr']=self.tau_fr    
        params['alpha']=self.alpha
        params['m_iso_mode']=self.m_iso_mode
        params['t_s']=self.t_s
        params['atmos']=self.atmos
        params['case']=self.case
        return json.dumps(params)
    
    #inverse of __str__
    def fromstr(self,_str):
        d = json.loads(_str)
        self.tau_fr=d['tau_fr']
        self.alpha=d['alpha']
        self.m_iso_mode=d['m_iso_mode']
        self.t_s=d['t_s']
        self.atmos=d['atmos']
        self.case=d['case']

    def __hash__(self):
        return hash((self.alpha,self.m_iso_mode,self.tau_fr,self.t_s,self.atmos,self.case))

    def __eq__(self,other):
        return (self.alpha,self.m_iso_mode,self.tau_fr,self.t_s,self.atmos,self.case) \
        == (other.alpha,other.m_iso_mode,other.tau_fr,other.t_s,other.atmos,other.case)


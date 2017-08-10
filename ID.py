import json

def sci(num):
    num = float(num)
    return "{:.1e}".format(num)

class Planet_ID: 

    def __init__(self,
                pebble_size : float = None,
                alpha : float = None,
                m_iso_mode=None,
                t_s=None,
                atmos:'clean or dusty'=None,
                init_a=None,
                case=None,
                tag=None,
                _dict:'pass in a dictionary containing all relevant parameter values for instant init'=None
                ):

        #used with optional dictionary normalization
        if _dict!=None:
            alpha= _dict['alpha']
            m_iso_mode= _dict['m_iso_mode']
            pebble_size= _dict['pebble_size']
            t_s= _dict['t_s']
            atmos=_dict['atmos']
            init_a=_dict['init_a']
            case=_dict['case']
            tag=_dict['tag']

        #below is the standard non-dict initialization
        self.alpha=sci(float(alpha))
        self.m_iso_mode=m_iso_mode
        self.pebble_size=sci(float(pebble_size))
        self.t_s=t_s
        self.atmos = atmos
        self.init_a = init_a
        self.case=case
        
        #rename this tag later for easier graphing
        self.tag=tag 

    #outputs the parameter values into a broken string that can be added to plots
    def param_list(self):

        ret =('alpha: '+ str(sci(self.alpha)) +
        '\npeb_size: '+str(sci(self.pebble_size)) +
        '\nt_s: '+str(self.t_s) +
        '\nm_iso_mode: '+str(self.m_iso_mode) +
        '\natmos: '+str(self.atmos) +
        '\ninit_a: '+str(self.init_a) +
        '\ncase: '+str(self.case))
           
        return ret

    #outputs the parameter values into a string that can be used for naming
    def file_name(self):

        return ('alpha '+str(sci(self.alpha))+
        ' peb_size '+str(sci(self.pebble_size))+
        ' t_s '+str(self.t_s)+
        ' m_iso_mode '+str(self.m_iso_mode)+
        ' atmos '+str(self.atmos)+
        ' init_a '+str(self.init_a)+
        ' case '+str(self.case))
        
    def __str__(self):
        params = {}
        params['pebble_size']=self.pebble_size    
        params['alpha']=self.alpha
        params['m_iso_mode']=self.m_iso_mode
        params['t_s']=self.t_s
        params['atmos']=self.atmos
        params['init_a']=self.init_a
        params['case']=self.case
        params['tag']=self.tag
        return json.dumps(params)
    
    #inverse of __str__
    def fromstr(self,_str):
        d = json.loads(_str)
        self.pebble_size=d['pebble_size']
        self.alpha=d['alpha']
        self.m_iso_mode=d['m_iso_mode']
        self.t_s=d['t_s']
        self.atmos=d['atmos']
        self.init_a=d['init_a']
        self.case=d['case']
        self.tag=d['tag']

    def __hash__(self):
        return hash((self.alpha,self.m_iso_mode,self.pebble_size,self.t_s,self.atmos,self.init_a,self.case,self.tag))

    def __eq__(self,other):
        return (self.alpha,self.m_iso_mode,self.pebble_size,self.t_s,self.atmos,self.init_a,self.case,self.tag) \
        == (other.alpha,other.m_iso_mode,other.pebble_size,other.t_s,other.atmos,other.init_a,other.case,other.tag)


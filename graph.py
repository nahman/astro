import numpy as np
import constants as ct
import disk as d 
import planet as pl
import ID
import json
import matplotlib.pyplot as plt


def load(file_name:str):
    loaded_dict={}
    with open(file_name) as json_data:
        d = json.load(json_data)
        for key, value in d.items():
  
            id_dict = json.loads(key)
            #id_dict['t_dep']=None
            newid=ID.Planet_ID(_dict=id_dict)
            loaded_dict[newid]=value
    return loaded_dict


#sol_dict is a keyed value in the dictionary outputted by solve/to_dict (keyed by planet id)
#sol_dict is a dict itself, keying time t by 't', core mass by 'm_core', atmosphere mass by 'm_atm', orbital distance by 'a'
#change var_of_interest to 'm_core', 'm_atm', or 'a' for a graph of only that variable

def graph(sol_dict,var_of_interest = 'all'): 
    colors = ['b','g','r','c','m','y']
    if var_of_interest!='all':
        plt.plot(sol_dict['t'], sol_dict[var_of_interest], label=var_of_interest)
        return 

    for key,value in sol_dict.items():
        if key!='t' and type(value)==list:
            plt.plot(sol_dict['t'],sol_dict[key],label=key)

def format_graph(title=None,xlabel=None,ylabel=None,log=True,t_stop=0):
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
        plt.title(title,fontsize=14)
    
        

    plt.tight_layout()
    #plt.show()

#in keeping with previous convention, the variable is what is plotted, the parameter is what is being changed. I assume in the sol_dict 
#only one param is being changed. Can't see as of now why you'd ever need to change two at a time.

#sols_dict a dictionary of many individual solutions. contrast w/ sol_dict, which is the dictionary of one solution.

#This funciton creates a single panel graph containing multiple plots of a single variable. Each plot is different from the other by the change of the same param.
def multi_graph(sols_dict,var_name,param_name):
    for key, value in sols_dict.items():
        plt.plot(value['t'],value[var_name],label=param_name+' '+str(key.__dict__[param_name]))
'''

sols_dict=load('tau_study.txt')
multi_graph(sols_dict,'a','tau_fr')
format_graph(title='effect of tau on orbital distance')

'''
def period(a):

    return  2*np.pi/(ct.G*ct.Msun/(a*ct.AU_to_cm)**3)**.5*(1/60) * (1/60) * (1/24)

#Gap opening criterion outlined in DM13
def gap_mass(a,alpha):

    disk_mach = a*ct.AU_to_cm/d.Disk.H(a)
    m_sh=(ct.Msun/ct.Mearth)/(disk_mach)**3
    return m_sh*17*(alpha*disk_mach)**.5

#plots mass of planet versus period
def p_plot(planet_id:ID.Planet_ID,sol_dict): 

    p_list = np.vectorize(period)(sol_dict['a'])
    plt.plot(p_list, [x+y for x,y in zip(sol_dict['m_core'],sol_dict['m_atm'])],label='m_p',color='k')
    plt.plot(p_list,np.vectorize(gap_mass)(sol_dict['a'],planet_id.alpha),label='gap-opening mass',linestyle='dashed',color='r')

    #point annotations
    #plt.axvline(x=period(sol_dict['a_stop']),linestyle='dashed',color='r',label='iso mass reached')
    #plt.plot(period(sol_dict['a_stop']),planet_id.m_iso_mode,'rx',label='m_iso')

#saves all graphs in the given file to the current directory
def save_graphs(filename:str):
    _d = load(filename)
    for key, value in _d.items():        
        f = plt.figure()
        p_plot(key,value)
        #graph(value)
        format_graph(title='planet mass vs. period',xlabel='orbital period, days',ylabel='planet mass, mearth')

        plt.figtext(1,.85,'params',fontsize=13)
        plt.figtext(1,.6,key.param_list())

        f.savefig(key.file_name()+'.pdf', bbox_inches='tight')
        plt.close()
        
    print('done!')

save_graphs('gap_mass_test.txt')


'''
y0=[1e-2,0,1]
t=np.logspace(-3,1,1000)



_d = load('case_2.txt')

for key,value in _d.items():
    #print(value['GCR'])
    plt.plot(value['t'], [x/y for x,y in zip(value['m_atm'],value['m_core'])],label='m_atm/m_core')
    plt.plot(value['t'],value['GCR'],linestyle='dashed',label='GCR eqn')
    #graph(value)
    format_graph(title='core and atm growth at constant 1AU',xlabel='time, myr',ylabel='m_atm/m_core')
    plt.show()


'''


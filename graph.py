import numpy as np
import constants as ct
import disk as d 
import planet as pl
import ID
import json
import matplotlib.pyplot as plt
from matplotlib import ticker

#same function as in solver_new
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


#sol_dict is a value in the dictionary outputted by solve/to_dict (keyed by planet id)
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

def format_graph(title=None,xlabel=None,ylabel=None,grid=False,log=True,t_stop=0):
    plt.legend(loc='best')
    
    if grid:
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

#period function used in period plotting
def period(a):
    return  2*np.pi/(ct.G*ct.Msun/(a*ct.AU_to_cm)**3)**.5*(1/60) * (1/60) * (1/24)

#Gap opening criterion outlined in DM13
def gap_mass(a,alpha):

    disk_mach = a*ct.AU_to_cm/d.Disk.H(a)
    m_sh=(ct.Msun/ct.Mearth)/(disk_mach)**3
    return m_sh*17*(alpha*disk_mach)**.5

#plots mass of planet versus period. gap mass overplotted.
def p_plot(planet_id:ID.Planet_ID,sol_dict:dict): 

    p_list = np.vectorize(period)(sol_dict['a'])
    plt.loglog(p_list, [x+y for x,y in zip(sol_dict['m_core'],sol_dict['m_atm'])],label='m_p',color='k')
    plt.loglog(p_list,np.vectorize(gap_mass)(sol_dict['a'],planet_id.alpha),label='gap-opening mass',linestyle='dashed',color='r')

#plots variables against orbital distance. vertical lines measure out the passage of time
def a_plot(planet_id:ID.Planet_ID,sol_dict:dict):

    index=0
    for time in sol_dict['t']:
        if index%10==0:
            plt.axvline(x=sol_dict['a'][index],linewidth=.5,color='.85')
        index+=1

    plt.loglog(sol_dict['a'], [x+y for x,y in zip(sol_dict['m_core'],sol_dict['m_atm'])],label='planet mass',color='.5')
    #plt.loglog(sol_dict['a'],np.vectorize(gap_mass)(sol_dict['a'],planet_id.alpha),label='gap-opening mass',color='r',linewidth=.5)
    plt.loglog(sol_dict['a'],sol_dict['m_core'],label='m_core',color='k',linestyle='dotted')
    plt.loglog(sol_dict['a'],sol_dict['m_atm'],label='m_atm',color='k',linestyle='dashed')
    plt.ylim(1e-2,1e+3)

#saves all graphs in the given file to the current directory
def save_graphs(filename:str,plot_mode: "set to p_plot to enforce period plotting"=None):
    _d = load(filename)
    for key, value in _d.items():        
        f,ax = plt.subplots()
        
        if plot_mode=='p_plot':
            p_plot(key,value)
            ax.xaxis.set_minor_formatter(ticker.NullFormatter())
            start, end = np.around(plt.xlim()[0],-1), np.round(plt.xlim()[1],-1)

            #manual x ticks. the stepsize can be changed
            plt.xticks(np.arange(start, end, 20))
            format_graph(title='planet mass vs period',xlabel='period, days',ylabel='planet mass, mearth')
        
        elif plot_mode=='a_plot':
            a_plot(key,value)

            ax.xaxis.set_minor_formatter(ticker.NullFormatter())
            start, end = np.around(plt.xlim()[0],1), np.around(plt.xlim()[1],1)

            #manual x ticks. the stepsize can be changed
            plt.xticks(np.arange(start, end,step=.05))
            format_graph(title='masses vs orbital dist',xlabel='dist, AU',ylabel='mass, mearth')

        else:
            graph(value) 
            format_graph(title='params vs time',xlabel='time, myr',ylabel='mass, mearth or dist, AU')

        plt.figtext(1,.85,'params',fontsize=13)
        plt.figtext(1,.6,key.param_list())

        f.savefig(key.file_name()+'.pdf', bbox_inches='tight')
        plt.close()
        
    print('done!')

#creates a three panel plot showing m_core, m_atm, a (one per panel). Each panel contains plots from multiple solutions, for comparison
def three_panel(filename:str,paramname:str,title:str):
    _d = load(filename)

    fig = plt.figure()

    m_core_ax = fig.add_subplot(3,1,1)
    m_atm_ax = fig.add_subplot(3,1,2,sharex=m_core_ax)
    a_ax = fig.add_subplot(3,1,3,sharex=m_core_ax)

    i=0
    colors = ['b','g','r','y','m','c']

    params_list=None

    for key, value in _d.items():
        params_list=key.param_list()
        col = colors[i%7]
        m_core_ax.loglog(value['t'],value['m_core'],col)
        m_atm_ax.loglog(value['t'],value['m_atm'],col)
        a_ax.loglog(value['t'],value['a'],col,label = paramname+': '+str(key.__dict__[paramname]))
        i+=1
    
    m_core_ax.set_title('m_core')
    m_core_ax.set_ylabel('mass, mearth')
    m_core_ax.grid()

    m_atm_ax.set_title('m_atm')
    m_atm_ax.set_ylabel('mass, mearth')
    m_atm_ax.grid()

    a_ax.set_title('orbital dist')
    a_ax.set_ylabel('dist, AU')
    a_ax.set_xlabel('time, Myr')
    a_ax.grid()

    a_ax.legend(loc='best')
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

    fig.suptitle(title)

    plt.text(18,2,params_list)
    #note that params_list needs to be edited to omit the paramter being varied
    
    plt.show()





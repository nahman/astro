import numpy as np
import constants as ct
import disk as d 
import planet as pl
import ID
import json
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.ticker import MaxNLocator

import csv

def sci(num):
    num = float(num)
    return "{:.1e}".format(num)

#same function as in solver_new
def load(file_name:str):
    loaded_dict={}
    t=None
    with open(file_name) as json_data:
        d = json.load(json_data)
        for key, value in d.items():
            if key=='t':
                t=value
                continue
            id_dict = json.loads(key)
            #id_dict['t_dep']=None
            newid=ID.Planet_ID(_dict=id_dict)
            loaded_dict[newid]=value
    return (loaded_dict,t)


#sol_dict is a value in the dictionary outputted by solve/to_dict (keyed by planet id)
#sol_dict is a dict itself, keying time t by 't', core mass by 'm_core', atmosphere mass by 'm_atm', orbital distance by 'a'
#change var_of_interest to 'm_core', 'm_atm', or 'a' for a graph of only that variable

def graph(sol_dict,var_of_interest = 'all',show_background=False): 
    if var_of_interest!='all':
        plt.plot(sol_dict['t'], sol_dict[var_of_interest], label=var_of_interest)
        return 
    
    t= sol_dict['t']
    plt.plot(t,sol_dict['m_core'],color='r',label='m_core')
    plt.plot(t,sol_dict['m_atm'],color='c',label='m_atm')
    plt.plot(t,sol_dict['a'],color='orange',label='a')

    if show_background:
        for key, val in sol_dict.items():
            if type(val) is list and key not in ['m_core','m_atm','a','t']:
                plt.plot(t,val,label=key)
    '''
    for key,value in sol_dict.items():
        if key!='t' and type(value)==list:
            plt.plot(sol_dict['t'],sol_dict[key],label=key)
    '''

def format_graph(title=None,xlabel=None,ylabel=None,grid=False,log=True,t_stop=0):
    plt.legend(bbox_to_anchor=(1.2, .4),
           bbox_transform=plt.gcf().transFigure)
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
    plt.ylim(ymin=1e-4)

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
    alpha = float(alpha)
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
    plt.loglog(sol_dict['a'],np.vectorize(gap_mass)(sol_dict['a'],planet_id.alpha),label='gap-opening mass',color='r',linewidth=.5)
    plt.loglog(sol_dict['a'],sol_dict['m_core'],label='m_core',color='k',linestyle='dotted')
    plt.loglog(sol_dict['a'],sol_dict['m_atm'],label='m_atm',color='k',linestyle='dashed')
    plt.ylim(1e-2,1e+3)

#saves all graphs in the given file to the current directory
def save_graphs(filename:str,plot_mode: "set to p_plot to enforce period plotting"=None,t=None):
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
            if value['died']:
                plt.axvline(x=value['t_death'],color='k',linestyle='dashed',label='t_death')
            format_graph(title='masses, dist vs time',xlabel='time, myr',ylabel='mass, mearth or dist, AU')
        '''
        plt.figtext(1,.85,'params',fontsize=13)
        if value['died']:
            plt.figtext(1,.6,key.param_list()+'\nplanet died')
        else:
            plt.figtext(1,.6,key.param_list())
        '''
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


def multpanel(key,sol_dict,t,save=False,loc_disk_m=False,show_dot=False):
    l_gap,l_jup,l_death = None, None, None

    fig = plt.figure(figsize=(8.5*.8,11*.8))
    fig.subplots_adjust(hspace=0)
    

    planet_ax = fig.add_subplot(5,1,1)
    a_ax = fig.add_subplot(5,1,2,sharex=planet_ax)
    disk_ax = fig.add_subplot(5,1,3,sharex=planet_ax)
    disk2_ax = fig.add_subplot(5,1,5,sharex=planet_ax)
    tau_ax = fig.add_subplot(5,1,4,sharex=planet_ax)

    if sol_dict['t_death']:
        l_death = disk2_ax.axvline(x=sol_dict['t_death'],color='k',linestyle='-.')
    if sol_dict['t_gap']:
        l_gap = disk2_ax.axvline(x=sol_dict['t_gap'],color='.75',linestyle='dashed')
    if sol_dict['t_jup']:
        l_jup = disk2_ax.axvline(x=sol_dict['t_jup'],color='.75',linestyle='dotted')


    if sol_dict['died']:
        a_ax.axvline(x=sol_dict['t_death'],color='k',linestyle='-.')
        planet_ax.axvline(x=sol_dict['t_death'],color='k',linestyle='-.')
        disk_ax.axvline(x=sol_dict['t_death'],color='k',linestyle='-.')
        tau_ax.axvline(x=sol_dict['t_death'],color='k',linestyle='-.')

    if sol_dict['t_gap']!=None:
        a_ax.axvline(x=sol_dict['t_gap'],color='.75',linestyle='dashed')
        planet_ax.axvline(x=sol_dict['t_gap'],color='.75',linestyle='dashed')
        disk_ax.axvline(x=sol_dict['t_gap'],color='.75',linestyle='dashed')
        tau_ax.axvline(x=sol_dict['t_gap'],color='.75',linestyle='dashed')
    
    if sol_dict['t_jup']!=None:
        a_ax.axvline(x=sol_dict['t_jup'],color='.75',linestyle='dotted')
        planet_ax.axvline(x=sol_dict['t_jup'],color='.75',linestyle='dotted')
        disk_ax.axvline(x=sol_dict['t_jup'],color='.75',linestyle='dotted')
        tau_ax.axvline(x=sol_dict['t_jup'],color='.75',linestyle='dotted')
    
    a_ax.axvline(x=key.t_s,color='k',linestyle='dotted')
    planet_ax.axvline(x=key.t_s,color='k',linestyle='dotted')
    disk_ax.axvline(x=key.t_s,color='k',linestyle='dotted')
    tau_ax.axvline(x=key.t_s, color='k',linestyle='dotted')
    l_ts = disk2_ax.axvline(x=key.t_s, color='k',linestyle='dotted')


    #if sol_dict['t_edge']!=None:
        #pass
        #a_ax.axvline(x=sol_dict['t_edge'],color='#013282',linestyle='dotted')
        #planet_ax.axvline(x=sol_dict['t_edge'],color='#013282',linestyle='dotted')
        #disk_ax.axvline(x=sol_dict['t_edge'],color='#013282',linestyle='dotted')
        #disk2_ax.axvline(x=sol_dict['t_edge'],color='#013282',linestyle='dotted',label='t_edge')
    
    #planet_ax.loglog(t,sol_dict['gap_mass'],color='.8',linestyle='dashed',label='gap_mass')
    #planet_ax.loglog(t,sol_dict['m_disk_s'],color='.5',linestyle='dotted',label='s_disk')
    planet_ax.loglog(t,sol_dict['m_disk'],color='#013282',linestyle='dotted',label=r'$M_g$')
    planet_ax.loglog(t,sol_dict['m_core'],color='.75',linestyle='-',label=r'$M_c$')
    planet_ax.loglog(t,sol_dict['m_atm'],color='.5',linestyle='-',label=r'$M_a$')

    planet_ax.set_ylim(ymin=1e-5)

    if not sol_dict['died']:
        final = sci(sol_dict['a'][-1])
        a_ax.axhline(y=sol_dict['a'][-1],color='.75',linestyle='dotted',label=r'$R_f$')
    
    a_ax.loglog(t,sol_dict['a'],color='k',label=r'$R$') 
    #a_ax.axhline(y=value['r1'],color='.25',linestyle='dotted',label=r'$R_1$')

    disk_ax.loglog(t,sol_dict['sigma_solid'],color='.5',linestyle='dotted',label=r'$\Sigma_{s}$')
    
    if sol_dict['t_gap']!=None:
        disk_ax.loglog(t,sol_dict['sigma_gas_p'],color='#013282',linestyle='-.',label=r'$\Sigma_{g}$')
    else:
        disk_ax.loglog(t,sol_dict['sigma_gas_up'],color='#013282',linestyle='-.',label=r'$\Sigma_{g}$')
    
    #disk_ax.loglog(t,sol_dict['mmen'],color='#013282',linestyle='dashed',label='sigma_mmen')

    tau_ax.loglog(t,sol_dict['tau'],label=r'$\tau_{fr}$',color='#994312',linestyle='dashed')
    tau_ax.axhline(y=float(key.alpha),label=r'$\alpha$',color='#994312',linestyle='dotted')

    disk2_ax.loglog(t,sol_dict['m_disk_dot'],label=r'$\dot M_{g}$',color='#013282',linestyle='dashed')

    planet_ax.legend(loc='best',frameon=False,framealpha=0,prop={'size':10})
    a_ax.legend(loc='best',frameon=False,framealpha=0,prop={'size':10})
    disk_ax.legend(loc='best',frameon=False,framealpha=0,prop={'size':10})
    tau_ax.legend(loc='best',frameon=False,framealpha=0,prop={'size':10})
    disk2_ax.legend(loc='best',frameon=False,framealpha=0,prop={'size':10})

    l = [[],[]]
    if l_gap != None:
        l[0].append(l_gap)
        l[1].append(r'$t_{gap}$')
    if l_jup != None:
        l[0].append(l_jup)
        l[1].append(r'$t_{jup}$')
    if l_death != None:
        l[0].append(l_death)
        l[1].append(r'$t_{end}$')
    l[0].append(l_ts)
    l[1].append(r'$t_s$')

    fig.legend(l[0],l[1],'lower right',prop={'size':10},frameon=False,framealpha=0)

    #planet_ax.text(1.01,.95,'params',transform=planet_ax.transAxes,fontsize=12)
    #planet_ax.text(1.01,-.6,key.param_list()+'\nm_core: '+str(sci(sol_dict['m_core'][-1]))+'\nR1: '+str(sci(sol_dict['r1'])),transform=planet_ax.transAxes)

    planet_ax.set_ylabel(r'$M/M_\bigoplus$')
    a_ax.set_ylabel(r'$R$/1 AU')
    disk_ax.set_ylabel(r'g/cm$^2$')
    disk2_ax.set_ylabel(r'$M_\bigodot$/yr')
    tau_ax.set_ylabel('(unitless)')

    
    tau_ax.tick_params(axis='x',which='both',bottom='off',labelbottom='off')

    disk2_ax.set_xlabel(r'$t$/1 Myr')
    
    if key.atmos=='clean':
        title = r'planet evolution for $s=$'+key.pebble_size+', dust-free envelope'
    else:
        title = r'planet evolution for $s=$'+key.pebble_size+', dusty envelope'
    plt.suptitle(title)

    if save:
        if key.tag!=None:
            fig.savefig(key.file_name()+' '+key.tag+'.pdf', bbox_inches='tight')
        else:
            fig.savefig(key.file_name()+'.pdf', bbox_inches='tight')
        plt.close()

_d,t=load('all.txt')
'''
for key, value in _d.items():
    if value['t_jup']==None:
        continue
    multpanel(key,value,value['t'],show_dot=True,save=True)
    break
plt.show()
'''

fig = plt.figure()
data=[]

for key, value in _d.items():
    if value['t_jup']!=None:
        data.append([period(value['a'][0]),period(value['a'][-1]),key.alpha,key.pebble_size])

for pi,pf,al,peb in data:
    if float(peb) == .1:
        color='yellow'
    elif float(peb) == 1:
        color='orange'
    else:
        color='red'
    if float(al) == 1e-2:
        shape='x'
    else:
        shape='.'
    plt.scatter(pi,pf,c=color,marker=shape)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$P_i$ / 1 day')
plt.ylabel(r'$P_f$ / 1 day')

plt.text(20,2000,r'$s$'+ '\n .1 \n 1 \n 10',multialignment='right')
plt.text(20,200,r'$\alpha$'+'\n 1e-2 \n 1e-3')

plt.title('Initial vs. Final Period of All Runs')


fig.savefig('scatter.pdf',figsize=(8,8))



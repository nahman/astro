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
            newid=ID.Planet_ID(_dict=id_dict)
            loaded_dict[newid]=value
    return loaded_dict


#sol_dict is a keyed value in the dictionary outputted by solve/to_dict (keyed by planet id)
#sol_dict is a dict itself, keying time t by 't', core mass by 'm_core', atmosphere mass by 'm_atm', orbital distance by 'a'
#change var_of_interest to 'mc', 'ma', or 'a' for a graph of only that variable
def graph(sol_dict,var_of_interest = 'all'): 
    colors = ['b','g','r','c','m','y']
    if var_of_interest!='all':
        plt.plot(sol_dict['t'], sol_dict[var_of_interest], label=var_of_interest)
        return 

    for key,value in sol_dict.items():
        if key!='t':
            plt.plot(sol_dict['t'],sol_dict[key],label=key)

def format_graph(title=None,xlabel=None,ylabel=None,showline=False,log=True,t_stop=0):
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
        plt.title(title)
    if showline:
        plt.axvline(x= t_stop) #marks when gap forms
    plt.tight_layout()
    plt.show()


d = load('boop.txt')
for key, value in d.items():
    if key.m_iso == 15 and value['a'][0]==1:
        graph(value)
        format_graph(title=key.name())


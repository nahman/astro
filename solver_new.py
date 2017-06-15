import numpy as np
import constants as ct
import disk as d 
import planet as pl
import ID

from scipy.integrate import odeint
import json

import os.path
from os import remove
import copy

toggled=False

#pass in a planet. returns diff system using said planet, for use in odeint
def diff_wrapper(_planet:pl.Planet):
    
    def diff_system(y,t):
        nonlocal _planet

        if _planet.a<=1e-2 or y[2]<=1e-2:
            global toggled
            if not toggled:
                print('planet died')
                toggled=True
            _planet.a=1e-2  #set a min value for orbital distance
            y[2]=1e-2
            
            return [0,0,0] 
        
        _planet.redefine(y[0],y[1],y[2],t)

        a_dot = _planet.a_dot(t)
        m_core_dot = _planet.m_core_dot(t)
        m_atm_dot = _planet.m_atm_dot(t)

        return [m_core_dot,m_atm_dot,a_dot]

    return diff_system
    
def solve(t,planet:pl.Planet,info=False):
    if info:
        sol,info = odeint(diff_wrapper(planet),planet.state(),t, full_output=True)
        return (t,sol,info)

    sol = odeint(diff_wrapper(planet),planet.state(),t)
    return (t,sol)

#check if the solution dictionary makes sense
def sol_check(sol:dict):
    min_c=1e-3
    max_c=1e3

    min_atm=0
    max_atm=1e3

    min_a=0
    max_a=20

    _l=zip(list(sol[:,0]),list(sol[:,1]),list(sol[:,2]))
    for x , y , z in _l:
        if not min_c<x<max_c or not min_atm<=y<max_atm or not min_a<z<max_a:
            print('possible error: '+'mc '+str(x)+' ma '+str(y)+' a '+str(z))
            return False
    return True
    
#Very experimental
#might not be catching all 'excess work done' cases
def wiggle_solve(t,planet:pl.Planet,max_iter=5):
    planetcopy=copy.deepcopy(planet)
    sol = solve(t,planetcopy)

    alpha=planet.disk.alpha
    alpha_exp = np.floor(np.log10(alpha))
    tau=planet.disk.tau_fr
    tau_exp= np.floor(np.log10(tau))

    attempt=0
    while not sol_check(sol[1]):
        if(attempt>max_iter):
            raise Exception('unable to solve')
        print('wiggling')
        planetcopy=copy.deepcopy(planet)
        planetcopy.disk.alpha=np.random.randint(-99,99)*(10**(alpha_exp-5))+alpha
        planetcopy.disk.tau_fr=np.random.randint(-99,99)*(10**(tau_exp-5))+tau
        
        #print('alpha: '+str(planetcopy.disk.alpha))
        #print('tau: '+str(planetcopy.disk.tau_fr))
        
        sol = solve(t,planetcopy)
        attempt+=1
    return sol

#Convert planet - solution pair to single entry dictionary keyed by planet_id. Makes writing to Json easier.
#Probably unecessary saving t. each file should really share the same t.

def to_dict(planet:pl.Planet,sol):
    vals={}
    vals['t']=list(sol[0])
    vals['m_core']=list(sol[1][:,0])
    vals['m_atm']=list(sol[1][:,1])
    vals['a']=list(sol[1][:,2])

    #NOT CGS
    vals['sigma_gas']= [x*ct.AU_to_cm**2/ct.Mearth for x in list(np.vectorize(planet.disk.new_Sigma_gas_t)(vals['a'],vals['t']))]

    #now for some points of interest
    vals['t_stop']=planet.tstop
    vals['a_stop']=planet.astop

    d={}
    d[planet.id]=vals
    return d

#save dictionary in json format
def write(d:dict,file_name:str):
    new_dict={} 
    with open(file_name, 'w') as fp:
        for key, value in d.items():
            new_dict[str(key)]=value
        json.dump(new_dict, fp,indent=4)
    return json.dumps(new_dict,indent=4)

#feed in planet, solution, filename. adds the new solution to the file, returns the new solution dictionary
#creates a new file if none exists
def addto(planet:pl.Planet,sol,file_name:str):

    new_dict=to_dict(planet,sol)
    if os.path.isfile(file_name):
        loaded_dict=load(file_name)

        loaded_dict.update(new_dict)
        write(loaded_dict,file_name)
        return loaded_dict
    write(new_dict,file_name)
    return new_dict
    
#loads given file returning the dictionary the file has stored.
def load(file_name:str):
    loaded_dict={}
    with open(file_name) as json_data:
        d = json.load(json_data)
        for key, value in d.items():
            id_dict = json.loads(key)
            newid=ID.Planet_ID(_dict=id_dict)
            loaded_dict[newid]=value
    return loaded_dict

def delete(filename:str):
    try:
        os.remove(filename)
    except:
        pass


#Case 1 testing
alpha_vals=[1.00013e-2,1.000012e-3,1.000012e-4,1.000012e-5]
tau_vals=[1.00015e-2,1.00012e-1,1.00011,1.0012e+1,1.00011e+2]
a_vals=[2,3,4,5,6,7,8,9,10]

y0=[1e-2,0,1]
t1=np.logspace(-3,1,1000)

#disk = d.Disk(1e-1,1e-4,1)
#planet = pl.Planet(y0[0],y0[1],y0[2],disk,case=1)
#sol=solve(t1,planet)
#addto(planet,sol,'error.txt')

'''
for tau in tau_vals:
    for alpha in alpha_vals:
        disk = d.Disk(tau,alpha,1e-1)
        planet = pl.Planet(y0[0],y0[1],y0[2],disk,case=1)
        print(planet.id.file_name())
        sol=solve(t1,planet)
        toggled=False
        addto(planet,sol,'t_s_1e-1.txt')
'''


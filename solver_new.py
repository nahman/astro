import numpy as np
import constants as ct
import disk as d 
import planet as pl
from scipy.integrate import odeint
import json
import ID



def diff_wrapper(planet:pl.Planet):
    
    def diff_system(y,t):
        nonlocal planet
        planet.redefine(y[0],y[1],y[2])
        m_core_dot = planet.m_core_dot(t)
        m_atm_dot = planet.m_atm_dot_dusty(t)
        a_dot = planet.a_dot(t)

        if planet.a<=1e-2:
            planet.a=1e-2
            a_dot=0

        return [m_core_dot,m_atm_dot,a_dot]

    return diff_system
    
def solve(t,planet:pl.Planet):
    sol = odeint(diff_wrapper(planet),planet.state(),t)
    return (t,sol)



#Convert planet - solution pair to single entry dictionary keyed by planet_id. Makes writing to Json easier.
#Probably unecessary saving t. each file should really share the same t.
def to_dict(planet:pl.Planet,sol):
    vals={}
    vals['t']=list(sol[0])
    vals['m_core']=list(sol[1][:,0])
    vals['m_atm']=list(sol[1][:,1])
    vals['a']=list(sol[1][:,2])
    d={}
    d[planet.id]=vals
    return d

#feed in planet, solution, filename. adds the new solution to the file, returns the new solution dictionary
def addto(planet:pl.Planet,sol,file_name:str):
    loaded_dict=load(file_name)
    new_dict = {}
    new_dict=to_dict(planet,sol)
    loaded_dict.update(new_dict)
    write(loaded_dict,file_name)
    return loaded_dict
    
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

#save dictionary in json format
def write(d:dict,name:str):
    new_dict={} 
    with open(name, 'w') as fp:
        for key, value in d.items():
            new_dict[str(key)]=value
        json.dump(new_dict, fp,indent=4)
    return json.dumps(new_dict,indent=4)




#example usage
'''
y0=[1e-2,0,1]
t=np.logspace(-3,1,100)

disk=d.Disk('pow',1,1e-4,t_init=1e-3)
planet =pl.Planet(y0[0],y0[1],y0[2],10,disk)

disk2=d.Disk('exp',1,1e-4,t_init=1e-3)
planet2=pl.Planet(y0[0],y0[1],y0[2],10,disk2)

sol =solve(t,planet)

write(to_dict(planet,sol),'boop.txt')

sol2 = solve(t,planet2)

addto(planet2,sol2,'boop.txt')
'''
t=np.logspace(-3,1,100)
disk=d.Disk('pow',1,1e-4,t_init=1e-3)
planet4 = pl.Planet(1e-3,0,2,15,disk)
sol4=solve(t,planet3)
addto(planet4,sol4,'boop.txt')

import numpy as np
import constants as ct
import disk as d 
import planet as pl
from scipy.integrate import odeint
import json
import ID
import os.path



def diff_wrapper(_planet:pl.Planet):
    
    def diff_system(y,t):
        nonlocal _planet

        _planet.redefine(y[0],y[1],y[2],t)

        if _planet.a<=1e-2:
            _planet.a=1e-2  #set a min value for orbital distance
            a_dot=0
        else:
            a_dot = _planet.a_dot(t)
        
        #print(_planet.a)

        m_core_dot = _planet.m_core_dot(t)
        m_atm_dot = _planet.m_atm_dot(t)

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

    #some extra info


    #now for some points of interest
    vals['t_stop']=planet.tstop
    vals['a_stop']=planet.astop

    d={}
    d[planet.id]=vals
    return d

#feed in planet, solution, filename. adds the new solution to the file, returns the new solution dictionary
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

#save dictionary in json format
def write(d:dict,file_name:str):
    new_dict={} 
    with open(file_name, 'w') as fp:
        for key, value in d.items():
            new_dict[str(key)]=value
        json.dump(new_dict, fp,indent=4)
    return json.dumps(new_dict,indent=4)



#M_Iso:5,10,15,20,25
#Tau:.1
#t_ini: -3,-2,-1
#dep_time: .1,1,10 #still altering this?
#alpha: -2,-3,-4,-5


'''
dep_time_vals=[.1,1,10]
m_iso_vals=[5.001,10.001,15.001,20.001] 
t_ini_vals=[1e-3,1e-2,1e-1]
alpha_vals=[1e-2,1e-3,1e-4,1e-5]

y0=[1e-2,0,1]
t1=np.logspace(-3,1,1000)
t2=np.logspace(-2,1,1000)
t3=np.logspace(-1,1,1000)
i=0
for m_iso in m_iso_vals:
    for alpha in alpha_vals:
        disk = d.Disk('mix',.1,alpha,t_init=1e-3)
        planet = pl.Planet(y0[0],y0[1],y0[2],m_iso,disk)
        print(planet.id.name())
        sol = solve(t1,planet)
        addto(planet,sol,'more.txt')
    i+=1
'''

'''

y0=[1e-2,0,1]
t=np.logspace(-3,1,1000)

disk = d.Disk('mix',.1,1e-5,1e-3)

planet = pl.Planet(y0[0],y0[1],y0[2],5,disk)


sol = solve(t,planet)

write(to_dict(planet,sol),'case_2.txt')
'''

y0=[1e-2,0,1]
t=np.logspace(-3,1,1000)

disk2 = d.Disk('mix',.1,1.001e-5,1e-3,atmos='clean')

planet2 = pl.Planet(y0[0],y0[1],y0[2],disk2)

sol2=solve(t,planet2)
addto(planet2,sol2,'gap_mass_test.txt')
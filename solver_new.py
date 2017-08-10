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

import matplotlib.pyplot as plt

toggled=False


#pass in a planet. returns diff system using said planet, for use in odeint
def diff_wrapper(_planet:pl.Planet):
    
    def diff_system(y,t):

        nonlocal _planet

        if y[2]<=1e-2:
            if not _planet.died:
                #print('planet died: '+str(t))
                _planet.died=True
                _planet.t_death=t
            _planet.a=1e-2  #set a min value for orbital distance
            y[2]=1e-2

            return [0,0,0] 

        elif _planet.died and _planet.a>1e-2:
            _planet.died=False
        #print(_planet.disk.disk_mass(t,y[2]),y[1],y[0])
        _planet.redefine(y[0],y[1],y[2],t)

        a_dot = _planet._a_dot
        m_core_dot = _planet._mc_dot
        m_atm_dot = _planet._ma_dot

        #stability enhancement
        if m_core_dot>1e10 or m_core_dot<0:
            m_core_dot=0
        if m_atm_dot>1e10 or m_atm_dot<0:
            m_atm_dot=0

        #print(m_atm_dot)
        #print(_planet.m_core,_planet.m_atm,_planet.a)
        return [m_core_dot,m_atm_dot,a_dot]



    return diff_system
    
def solve(t,planet:pl.Planet,info=False):

    if planet.disk.R1<planet.a:
        print('ignored')
        return

    if info:
        sol,info = odeint(diff_wrapper(planet),planet.state(),t, full_output=True)
        return (t,sol,info)

    sol = odeint(diff_wrapper(planet),planet.state(),t)
    return (t,sol)

#check if the solution dictionary makes sense
def sol_check(sol:dict):
    min_c=.8 * sol[:,0][0]
    max_c=5e2

    min_atm=-1e-6
    max_atm=4e2

    min_a=0
    max_a=11

    _l=zip(list(sol[:,0]),list(sol[:,1]),list(sol[:,2]))
    for x , y , z in _l:
        if not min_c<x<max_c or not min_atm<= y<max_atm or not min_a<z<max_a:
            print('possible error: '+'mc '+str(x)+' ma '+str(y)+' a '+str(z))
            return False
    return True
    
#Very experimental
#Seems to be rather robust though

#solves the system, then checks if the output makes sense. if not, retry solving the system,
#after altering alpha and pebble size by insignificant amounts
def wiggle_solve(t,planet:pl.Planet,max_iter=5,accuracy=5):
    planetcopy=copy.deepcopy(planet)
    sol = solve(t,planetcopy)
    if sol == None:
        return
    alpha=planet.disk.alpha
    alpha_exp = np.floor(np.log10(alpha))
    peb=planet.disk.pebble_size
    peb_exp=np.floor(np.log10(peb))

    attempt=0
    while not sol_check(sol[1]):
        if(attempt>max_iter):
            raise Exception('unable to solve')
        print('planet wiggling...')
        planetcopy=copy.deepcopy(planet)
        planetcopy.disk.alpha=np.random.randint(-99,99)*(10**(alpha_exp-accuracy))+alpha
        planetcopy.disk.pebble_size=np.random.randint(-99,99)*(10**(peb_exp-accuracy))+peb
        
        sol = solve(t,planetcopy)
        attempt+=1
    
    #transfer parameters
    planet.died=planetcopy.died
    planet.t_death=planetcopy.t_death
    planet.t_stop=planetcopy.t_stop
    planet.t_gap=planetcopy.t_gap
    planet.t_jup=planetcopy.t_jup
    planet.disk.t_edge=planetcopy.disk.t_edge
    return sol

#Convert planet - solution pair to single entry dictionary keyed by planet_id. Makes writing to Json easier. Has lots of additional data.

'''
key                 value

m_core              core mass
m_atm               atmosphere mass
a                   orbital radius
sigma_gas_up        unperturbed gas surface density
sigma_gas_p         perturbed gas surface density (takes into account planet presence)
sigma_solid         surfaced density of dust/pebbles for core accretion
mmen                minimum mass extra solar nebula (gas)
m_disk              total mass of gas in disk outside of planet (current outer integration limit 100AU)
m_disk_s            same as m_disk but for solids
m_disk_dot          gas accretion rate by star of disk
tau                 dimensionless friction time of pebbles
gap_mass            planet mass necessary to open a disk of certain viscosity at certain location, etc
loc_disk_m          local gas disk mass

t_gap               time of gap opening in disk
died                boolean, whether the planet fell into the star or not
t_jup               time of runaway gas accretion and gas giant formation
t_death             time of planet death (via infall into the star)
t_edge              time at which outer edge of the solid disk sweeps past the planet as dust/pebbles drift inwards
r1                  dimensionless radial scale factor of the disk, for more info look in Hartmann 1998 <or maybe 1997?
'''

def to_dict(_planet:pl.Planet,sol):
    vals={}

    t_arr = list(sol[0])
    a_arr = list(sol[1][:,2])
    mc_arr = list(sol[1][:,0])
    ma_arr = list(sol[1][:,1])
    mp_arr = [x+y for x,y in zip(mc_arr,ma_arr)]

    #saving the main paramters
    vals['m_core']= mc_arr
    vals['m_atm']= ma_arr
    vals['a'] = a_arr
    vals['t'] = t_arr

    #surface density panel
    vals['sigma_gas_up'] = np.vectorize(_planet.disk.new_Sigma_gas_t)(a_arr,t_arr).tolist()

    if _planet.t_gap!=None:
        f = lambda a,t: _planet.disk.new_Sigma_gas_t(a,t) if t<=_planet.t_gap else _planet.disk.new_Sigma_gas_t(a,t,gap=True)
        vals['sigma_gas_p'] = np.vectorize(f)(a_arr,t_arr).tolist()
    
    vals['sigma_solid'] = np.vectorize(_planet.disk.sigma_solid)(a_arr,t_arr).tolist()
    vals['mmen'] = np.vectorize(_planet.disk.Sigma_mmen)(a_arr).tolist()

    #disk masses and such
    f = lambda t,a: _planet.disk.disk_mass(t,a)
    vals['m_disk'] = np.vectorize(f)(t_arr,a_arr).tolist()

    #f = lambda a, t: _planet.disk.disk_mass(t,a,mode='solid')
    #vals['m_disk_s'] = np.vectorize(f)(a_arr,t_arr).tolist()

    f = lambda a, t: _planet.disk.new_Sigma_gas_t(a,t)*(a*ct.AU_to_cm)**2/(ct.Mearth*ct.m_jup)
    vals['loc_disk_m'] = list(np.vectorize(f)(a_arr,t_arr))

    vals['m_disk_dot'] = np.vectorize(_planet.disk.disk_mass_dot)(t_arr).tolist()

    #misc
    vals['tau'] = np.vectorize(_planet.disk.tau)(a_arr,t_arr).tolist()
    vals['gap_mass'] = np.vectorize(_planet.gap_mass_alt)(a_arr).tolist()


    #extra points of interest
    vals['t_gap']=_planet.t_gap
    vals['died']=_planet.died
    vals['t_death']=_planet.t_death
    vals['r1']=_planet.disk.R1
    vals['t_jup']=_planet.t_jup
    
    if 0 in vals['sigma_solid']:
        vals['t_edge']=t_arr[vals['sigma_solid'].index(0)]
    else:
        vals['t_edge']=None

    d={}
    d[_planet.id]=vals
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
    if t==None:
        return loaded_dict
    return (loaded_dict,t)

def delete(filename:str):
    try:
        os.remove(filename)
    except:
        pass


_d1=load('rev_03.txt')
_d2=load('rev_1.txt')
_d3=load('rev_3.txt')
_d4=load('rev_5.txt')
_d5=load('rev_10.txt')

_d1.update(_d2)
_d1.update(_d3)
_d1.update(_d4)
_d1.update(_d5)


write(_d1,'all.txt')

'''
#edit these if you don't have a params.py
t_s = [1] 
alphas = [1e-2]
init_a = [3] #[.3,1,3,5,10]
atmos = ['clean']
peb_size = [10]
start = -3
end = 2
steps = 3500
cm = 1e-4
smass = 30

try:
    from params import *
except:
    print('params.py does not exist or has an error. using default values')


t=np.logspace(start,end,steps)
_dict={}

for peb in peb_size:
    for atm in atmos:
        for a in init_a:
            for ts in t_s:
                for alpha in alphas:
                    disk=d.disk(peb,alpha,ts,atmos=atm,s_mass=smass)
                    planet=pl.Planet(cm,0,a,disk)
                    print(planet.id.file_name())
                    sol=wiggle_solve(t,planet)
                    print('solved')
                    if sol !=None:
                        _dict.update(to_dict(planet,sol))


#_dict['t']=t.tolist()

write(_dict,'3error.txt')
'''
# this program is adopted from csvoutput.py and modified
# for binary file output instead of csv output.
# 12/4/2018

########### constants ##########
# constants are shared with n_body_simulation.py
# (do not change)
parsec = 3.086 * 10**16 #[m]
Msun = 1.989 * 10**30 #[kg]
G = 6.67408 * 10**(-11) #[m^3 * kg^(-1) * s^(-2)]
pi = 3.141592653
day = 3600*24 #[sec]
year = 3600*24*365 #[sec]
kpc = 3.086 * 10**19 #[m]

########### parameters ##########
# Data file
filename = 'initial_data/initial_data.dat'

# Galaxy property
gal_disk_r  = 25 * 10**3 * parsec #[m]
gal_disk_dz = 0.15 * 10**3 * parsec #[m]
gal_bulge_r = 0.5 * 10**3 * parsec #[m]
BH_mass = 8.2 * 10**36 #[kg]
rho_0 = 6.0 * 10 ** (-21) #[kg/m^3]
r_c = 60 * kpc #[m]

# Star property
star_v = 200 * 10**3 #[m/s]
num_stars   = 1* 10**2
actual_num  = 10 ** 10
mass_coef   = actual_num / num_stars

# simulation parameters
dt = 10**6 * year #[sec]
t_max = 10**9 * year #[sec]
####################################


def initial_list_generator():
    import numpy as np
    # position generation
    num = np.arange(num_stars)
    mass = np.random.uniform(1*Msun*mass_coef,20*Msun*mass_coef,num_stars)
    x = np.random.triangular(-1*gal_disk_r,0,gal_disk_r,num_stars)
    y = []
    for i in range(len(x)):
        ylim = np.sqrt(gal_disk_r**2 - x[i]**2)
        y.append(*np.random.uniform(-1*ylim,ylim,1))
    z = np.random.uniform(-1*gal_disk_dz,gal_disk_dz,num_stars)

    # velocity generation
    vx = []
    vy = []
    vz = []
    for i in range(num_stars):
        d = np.sqrt(x[i]**2 + y[i]**2)
        v = np.random.uniform(star_v - 150*10**3, star_v + 250*10**3,1)
        vx.append(float(star_v * -1 * y[i] / d))
        vy.append(float(star_v * x[i] / d))
        vz.append(float(np.random.uniform(-10*10**3,-10*10**3,1)))

    # starlist
    initial_list = []
    for i in range(num_stars):
        stardata = [num[i],mass[i],x[i],y[i],z[i],vx[i],vy[i],vz[i]]
        initial_list.append(stardata)
    return initial_list


def condition_data_generator():
    names  = ['disk_r','disk_dz','bulge_r','BH_m','DM_rho0','r_c','star_v',\
            'num_stars','actual_num','mass_coef','dt','t_max']
    values = [gal_disk_r,gal_disk_dz,gal_bulge_r,BH_mass,rho_0,r_c,star_v,\
            num_stars,actual_num,mass_coef,dt,t_max]
    condition_data = [names,values]
    print(condition_data)
    return condition_data
  
def output(filename,condition_data,initial_list):
    import pickle
    output = open(filename,'wb')
    outputdata = [condition_data,initial_list]
    pickle.dump(outputdata,output)
    output.close()
    return True


### main ###
initial_list   = initial_list_generator()
condition_data = condition_data_generator()
if output(filename,condition_data,initial_list):
    print('Data was successfully created and saved to {}'.format(filename))


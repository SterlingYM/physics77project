# this program is adopted from csvoutput.py and modified
# for binary file output instead of csv output.
# 12/4/2018

########### constants ##########
# constants are shared with n_body_simulation.py
# (do not change)
parsec      = 3.086 * 10**16 #[m]
Msun        = 1.989 * 10**30 #[kg]
G           = 6.67408 * 10**(-11) #[m^3 * kg^(-1) * s^(-2)]
pi          = 3.141592653
day         = 3600*24 #[sec]
year        = 3600*24*365 #[sec]
kpc         = 3.086 * 10**19 #[m]

########### parameters ##########
# Data file
filename = 'initial_data.dat'

# Galaxy property
gal_disk_r  = 25 * kpc #[m]
gal_disk_dz = 0.15 * kpc #[m]
gal_bulge_r = 10 * kpc #[m]
gass_mass   = 10**9 * Msun #[kg]
BH_mass     = 8.2 * 10**36 + gass_mass #[kg]
rho_0       = 4 * 10**7 * Msun / (kpc**3) #[kg/m^3] #from arXiv.1304.5127
r_c         = 60 * kpc #[m] #from arXiv:astro-ph/0403206

# Star property
star_v      = 150 * 10**3 #[m/s]
totalstar   = 300
bulge_coef  = 0.4
bulgestar   = int(float(bulge_coef * totalstar))
num_stars   = totalstar - bulgestar
actual_num  = 10 ** 11
mass_coef   = actual_num / num_stars

# simulation parameters
dt          = 1 * 10**5 * year #[sec]
t_max       = 5* 10**9 * year #[sec]
softening   = 10 * gal_disk_r / (num_stars)**(1/2) #[m]: mean distance
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
        v = np.random.uniform(star_v * 0.8, star_v * 1.2,1)
        vx.append(float(v * -1 * y[i] / d))
        vy.append(float(v * x[i] / d))
        vz.append(float(np.random.uniform(-10*10**3,10*10**3,1)))

    # starlist
    initial_list = []
    for i in range(num_stars):
        stardata = [num[i],mass[i],x[i],y[i],z[i],vx[i],vy[i],vz[i]]
        initial_list.append(stardata)
    return initial_list


def bulge ():
    import numpy as np
    # spherical coordinate
    num_bulge = int(bulgestar / 2)
    radius = np.random.uniform(0.0,gal_bulge_r,num_bulge)  
    theta = np.random.uniform(0.,1.,num_bulge)*pi
    phi = np.arccos(1-2*np.random.uniform(0.0,1.,num_bulge))
    
    # number index
    num = np.arange(num_bulge * 2) 

    # position
    x1 = radius * np.sin( theta ) * np.cos( phi )
    x  = np.concatenate((x1,-x1),axis=0)
    y1 = radius * np.sin( theta ) * np.sin( phi )
    y  = np.concatenate((y1,-y1),axis=0)
    z1 = radius * np.cos( theta )
    z  = np.concatenate((z1,z1),axis=0)

    # mass
    massbulge = np.random.uniform(1*Msun*mass_coef,20*Msun*mass_coef,bulgestar)

    # velocity
    vx = []
    vy = []
    vz = []
    for i in range(num_bulge*2):
        d = np.sqrt(x[i]**2 + y[i]**2)
        vel_net = np.sqrt(massbulge[i]*G/d)*(gal_bulge_r - abs(z[i]))/gal_bulge_r * 0.8 
        vx.append(float(vel_net * -1 * y[i] / d))
        vy.append(float(vel_net * x[i] / d))
        vz.append(float(np.random.uniform(-10*10**3,10*10**3,1)))
    
    # starlist_bulge
    starlist_bulge = []
    for i in range(num_bulge * 2):
        each_star = [num[i],massbulge[i],x[i],y[i],z[i],vx[i],vy[i],vz[i]]
        starlist_bulge.append(each_star)
    return starlist_bulge
    
    


def condition_data_generator():
    names  = ['disk_r','disk_dz','bulge_r','BH_m','DM_rho0','r_c','star_v',\
            'num_stars','bulgestars','actual_num','mass_coef','dt','t_max','softening']
    values = [gal_disk_r,gal_disk_dz,gal_bulge_r,BH_mass,rho_0,r_c,star_v,\
            num_stars,bulgestar,actual_num,mass_coef,dt,t_max,softening]
    condition_data = [names,values]
    print("Initial Condition Parameters:")
    for i in range(len(names)):
        print("\t* {:10}\t= {:1.2}".format(names[i],float(values[i])))

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
bulge_list     = bulge()
initial_list = initial_list + bulge_list
if output(filename,condition_data,initial_list):
    print('Data was successfully created and saved to \'{}\''.format(filename))


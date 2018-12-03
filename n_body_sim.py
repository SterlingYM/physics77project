#### constants ####
G = 6.67408 * 10**(-11) #[m^3 * kg^(-1) * s^(-2)]

#### parameters ####
# choose proper dt and t_max
csv_filename = 'Initial_Conditions_test.csv'
savefile_name = 'gal_hist.dat'
pi = 3.141592653
day = 3600*24 #[sec]
year = 3600*24*365 #[sec]
dt = 10**6 * year #[sec]
t_max = 10**9 * year #[sec]
BH_mass = 8.2 * 10**36 #[kg]
DM_density = 6.0 * 10 ** (-22) #[kg/m^3]

#### functions ####
def data_read(cev_filename):
    # imports csv file and returns initial condition data
    # formatted as [[star#,mass,x,y,vx,vy],[],[],...]
    initial_list = []
    import csv
    csvfile =  open(cev_filename, 'r')
    csvreader = csv.reader(csvfile, delimiter=',')
    for row in csvreader:
        initial_list.append(row)
    for i in range(len(initial_list)):
        initial_list[i][1] = float(initial_list[i][1])
        initial_list[i][2] = float(initial_list[i][2])
        initial_list[i][3] = float(initial_list[i][3])
        initial_list[i][4] = float(initial_list[i][4])
        initial_list[i][5] = float(initial_list[i][5])
        initial_list[i][6] = float(initial_list[i][6])
        initial_list[i][7] = float(initial_list[i][7])
        print('\rReading initial list...  \t{}/{} done'.format(i+1,len(initial_list)),end="")
    print("")
    return  initial_list

def time_development(initial_list,dt,t_max):
    # calls starloop() and stacks result on 'gal_hist'
    # repeats above for given time length
    gal_hist = []
    t = 0

    ## animation only (read animation block below) #####
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize = plt.figaspect(0.7)*1.5)
    ax = fig.add_subplot(111, projection='3d')
    ####################################################
    
    for k in range(int(t_max/dt)):
        if k == 0:
            starlist = starloop(initial_list)
            gal_hist.append([t,starlist])
        else:
            gal_hist.append([t,starloop(gal_hist[k-1][1])])

        ###  Real-time animation   #####################
        # To perform/disable real-time results, uncomment/comment below:
        animate_inline(gal_hist[k][1],ax,t) 
        ################################################

        t += dt
        print('\rSimulation in progress...\t{:.2f}% done'.format(float((k+1)*100)/float(t_max/dt)),end="")

    print("")
    return gal_hist

def starloop(starlist):
    # repeats calculation of new starlist for star[i]
    # calls net_force(), accel(), pos(), vel() within a loop for i
    import numpy as np
    
    new_starlist = []
    n_stars = len(starlist)
    
    # previous data
    p_pos = np.array([starlist[i][2:5] for i in range(n_stars)])
    p_vel = np.array([starlist[i][5:8] for i in range(n_stars)])
    # acceleration from previous data
    net_F = net_force(starlist)
    mass  = np.array([[starlist[i][1]] for i in range(n_stars)])
    acc   = np.divide(net_F,mass)

    # new data calculation
    pos   = p_pos + np.multiply(p_vel,dt) + np.multiply(acc,(dt**2)/2)
    vel   = p_vel + np.multiply(acc,dt)

    for i in range(len(starlist)):
        new_data = [i,*mass[i],*pos[i,:],*vel[i,:]]
        new_starlist.append(new_data)
    return new_starlist


def animate(gal_hist):
    # animates data passed in gal_hist
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize = plt.figaspect(0.7)*1.5)
    ax = fig.add_subplot(111, projection='3d')

    # boundary
    MAX = 0.8 * 10**21
    ax.set_xlim(-MAX,MAX)
    ax.set_ylim(-MAX,MAX)
    ax.set_zlim(-MAX,MAX)
    
    for i in range(len(gal_hist)):
        ax.cla()
        # data
        xlist = [gal_hist[i][1][j][2] for j in range(len(gal_hist[-1][1]))]
        ylist = [gal_hist[i][1][j][3] for j in range(len(gal_hist[-1][1]))]
        zlist = [gal_hist[i][1][j][4] for j in range(len(gal_hist[-1][1]))]
        
        #plot
        ax.scatter(0,0,0,color='orange')
        ax.scatter(xlist,ylist,zlist)
        plt.title('galaxy age = {:1.0} [yr]'.format(gal_hist[i][0]/year))
        plt.pause(0.001)
    plt.show()

def animate_inline(starlist,ax,t):
    # animates data real time
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    ax.cla()
    xlist = [starlist[i][2] for i in range(len(starlist))]
    ylist = [starlist[i][3] for i in range(len(starlist))]
    zlist = [starlist[i][4] for i in range(len(starlist))]
    ax.scatter(0,0,0,color='orange')
    ax.scatter(xlist,ylist,zlist)
    
    MAX = 0.8 * 10**21
    ax.set_xlim(-MAX,MAX)
    ax.set_ylim(-MAX,MAX)
    ax.set_zlim(-MAX,MAX)
 
    plt.title('galaxy age = {:1.3} [yr]'.format(t/year))
    plt.pause(0.001)

def dist(x1,x2,y1,y2,z1,z2):
    # calculates the distance between two sets of coordinates 
    import numpy as np
    dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return dist


def force_ij(star1,star2):
    # calculates the value of force between two particles given
    m1, x1, y1 ,z1 = star1[1], star1[2], star1[3], star1[4]
    m2, x2, y2, z2 = star2[1], star2[2], star2[3], star2[4]
    if x1==x2 and y1==y2 and z1==z2:
        return 0
    d = dist(x1,x2,y1,y2,z1,z2)
    force_ij = (G*m1*m2)/(d**2)
    Fx,Fy,Fz = components(force_ij,x1,x2,y1,y2,z1,z2)
    return [Fx,Fy,Fz]


def components(val,x1,x2,y1,y2,z1,z2):
    # returns given value in component form
    # assuming the value passed is a vector on (x1,y1) to (x2,y2)
    d = dist(x1,x2,y1,y2,z1,z2)
    val_x = val * ((x2 - x1)/d)
    val_y = val * ((y2 - y1)/d)
    val_z = val * ((z2 - z1)/d)
    return val_x, val_y, val_z


def net_force(starlist):
    # calls force_ij and takes the sum of all force in x & y direction
    # resulting value after for loop should be the net force on star[i]
    net_Fx, net_Fy, net_Fz = 0, 0, 0
    n_stars = len(starlist)

    import numpy as np
 
    # use multiprocessing and get all combination of forces
    
    import multiprocessing
    import itertools
    from itertools import product

    #with multiprocessing.Pool(n_stars) as pool:
    #    F_ij_combinations_xyz = pool.starmap(force_ij, product(starlist,repeat=2))
    arg1, arg2 = zip(*product(starlist,repeat=2))
    F_ij_combinations_xyz = list(map(force_ij,arg1,arg2))
    
    # change combination list to matrix
    F_ij_matrix = np.zeros((n_stars,n_stars,3),dtype=float)
    counter = 0
    for i in range(n_stars-1):
        for j in range(i+1,n_stars):
            F_ij_matrix[i,j,] = np.transpose(F_ij_combinations_xyz[counter]) # xyz-components
            counter += 1
    F_ij_transp = (-1)*F_ij_matrix.transpose((1,0,2))

    # sum all force on each object to get array of net force vector [[F_net_x],[F_net_y],[F_net_z]]
    net_F = np.sum((F_ij_matrix + F_ij_transp),axis=1)
    
    # add SMB gravity
    BH_data  = ['BH',BH_mass,0,0,0,0,0,0]
    BH_force = np.array([force_ij(starlist[i],BH_data) for i in range(n_stars)])
    
    # add DM gravity
    DM_force = np.array([force_ij(starlist[i],['DM',DM_mass(starlist[i]),0,0,0,0,0,0]) for i in range(n_stars)])
    
    return net_F + BH_force + DM_force

def DM_mass(star):
    r = dist(star[2],0,star[3],0,star[4],0)
    v = 4 * pi * r**3 / 3
    return v * DM_density

def save_data(gal_hist,savefile_name):
    # saves all data into .dat file
    import pickle
    output = open(savefile_name,'wb')
    pickle.dump(gal_hist,output)
    output.close()
    print('Generated Data was saved to {} successfully.'.format(savefile_name))

#### main ####
initial_list = data_read(csv_filename)
gal_hist = time_development(initial_list,dt,t_max)
save_data(gal_hist,savefile_name)

#animate(gal_hist) -> 12/2/2018: additional .py code will do the animation part now



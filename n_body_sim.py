#### constants ####
G = 6.67408 * 10**(-11) #[m^3 * kg^(-1) * s^(-2)]

#### parameters ####
# choose proper dt and t_max
csv_filename = 'Initial_Conditions_test.csv'
day = 3600*24 #[sec]
year = 3600*24*365 #[sec]
dt = 10000 * year #[sec]
t_max = 10**7 * year #[sec]
BH_mass = 8.2 * 10**36 #[kg]

#### functions ####
def data_read(cev_filename):
    # imports csv file and returns initial condition data
    # formatted as [[star#,mass,x,y,vx,vy],[],[],...]
    initial_list = []
    import csv
    csvfile =  open(cev_filename, 'r')
    csvreader = csv.reader(csvfile, delimiter=',')
    print('reading initial list... '),
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
        print('\r{}/{} done'.format(i+1,len(initial_list)),end="")
    print("")
    return  initial_list


def dist(x1,x2,y1,y2,z1,z2):
    # calculates the distance between two sets of coordinates 
    import numpy as np
    dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return dist


def force_ij(starlist1,starlist2):
    # calculates the value of force between two particles given
    m1, x1, y1 ,z1 = starlist1[1], starlist1[2], starlist1[3], starlist1[4]
    m2, x2, y2, z2 = starlist2[1], starlist2[2], starlist2[3], starlist2[4]
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

    # use multiprocessing and get all combination of forces
    import multiprocessing
    from itertools import product
    import numpy as np
    with multiprocessing.Pool(2) as pool:
        F_ij_combinations_xyz = pool.starmap(force_ij, product(starlist,repeat=2))
        # output will be [F_ij_x,F_ij_y,F_ij_z]
    
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

    #BH_Fx, BH_Fy, BH_Fz = components(starlist[i],['BH',BH_mass,0,0,0,0,0,0])
    #net_Fx += BH_Fx
    #net_Fy += BH_Fy
    #net_Fz += BH_Fz
    return net_F

def accel(mass,Fx,Fy,Fz):
    # calculates acceleration from net force ans mass
    # Newton's second law
    acc_x = Fx / mass
    acc_y = Fy / mass
    acc_z = Fz / mass
    return acc_x, acc_y, acc_z


def pos(starlist,i,ax,ay,az,dt):
    # calculates new position from current position, velocity, and acceleration
    _,_,x,y,z,vx,vy,vz = starlist[i]
    pos_x = x + vx*dt + ax*(dt**2)/2
    pos_y = y + vy*dt + ay*(dt**2)/2
    pos_z = z + vz*dt + az*(dt**2)/2
    return pos_x, pos_y, pos_z


def vel(starlist,i,ax,ay,az,dt):
    # calculates new velocity from current velocity and acceleration
    vx,vy,vz = starlist[i][5], starlist[i][6], starlist[i][7]
    v_x = vx + ax*dt
    v_y = vy + ay*dt
    v_z = vz + az*dt
    return v_x, v_y, v_z


def starloop(starlist):
    # repeats calculation of new starlist for star[i]
    # calls net_force(), accel(), pos(), vel() within a loop for i
    new_starlist = []
    n_stars = len(starlist)


    net_F = net_force(starlist)
    mass  = np.transpose([starlist[i][1] for i in range(n_stars)])
    acc   = np.divide(net_F,mass)
    pos   = 
    for i in range(len(starlist)):
        mass = starlist[i][1]
        Fx,Fy,Fz = net_F[i,0],net_F[i,1],net_F[i,2]
        ax,ay,az = accel(mass,Fx,Fy,Fz)
        x_new, y_new, z_new   = pos(starlist,i,ax,ay,az,dt)
        vx_new, vy_new, vz_new = vel(starlist,i,ax,ay,az,dt)
        new_starlist.append([starlist[i][0],mass,x_new,y_new,z_new,vx_new,vy_new,vz_new])
    return new_starlist

def time_development(initial_list,dt,t_max):
    # calls starloop() and stacks result on 'gal_hist'
    # repeats above for given time length
    gal_hist = []
    t = 0
    print('Simulation in progress... ')
    for k in range(int(t_max/dt)):
        if k == 0:
            starlist = starloop(initial_list)
            gal_hist.append([t,starlist])
        else:
            gal_hist.append([t,starloop(gal_hist[k-1][1])])
        t += dt
        #print(*gal_hist[k],sep="\n")
        print('\r{:.2f}% done'.format((k+1)*100/int(t_max/dt)),end="")
    print("")
    return gal_hist



#def animate():
    #animation part here 

#### main ####
initial_list = data_read(csv_filename)
gal_hist = time_development(initial_list,dt,t_max)
#animate(gal_hist)

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
fig = plt.figure(figsize = (15,15))
ax = fig.add_subplot(111, projection='3d')


for i in range(len(gal_hist)):
    ax.cla()
    #ax.set_xlim(-100,100)
    #ax.set_ylim(-100,100)
    #ax.set_zlim(-100,100)
    xlist = [gal_hist[i][1][j][2] for j in range(len(gal_hist[-1][1]))]
    ylist = [gal_hist[i][1][j][3] for j in range(len(gal_hist[-1][1]))]
    zlist = [gal_hist[i][1][j][4] for j in range(len(gal_hist[-1][1]))]
    ax.scatter(0,0,0,color='orange')
    ax.scatter(xlist,ylist,zlist)
    plt.title('galaxy age = {:.3f} [yr]'.format(gal_hist[i][0]/year))
    plt.pause(0.01)
    
plt.show()

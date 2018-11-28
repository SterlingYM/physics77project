#### constants ####
G = 6.67408 * 10**(-11) #[m^3 * kg^(-1) * s^(-2)]

#### parameters ####
# choose proper dt and t_max
csv_filename = 'InitialCondition.csv'
dt = 500 #[s]
t_max = 10**4

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

    return  initial_list


def dist(x1,x2,y1,y2):
    # calculates the distance between two sets of coordinates 
    import math
    dist = math.hypot(x2 - x1, y2 - y1)
    return dist


def force_ij(m1,m2,dist):
    # calculates the value of force between two particles given
    force_ij = (G*m1*m2)/(dist**2)
    return force_ij


def components(val,x1,x2,y1,y2):
    # returns given value in component form
    # assuming the value passed is a vector on (x1,y1) to (x2,y2)
    d = dist(x1,x2,y1,y2)
    val_x = val * ((x2 - x1)/d)
    val_y = val * ((y2 - y1)/d)
    return val_x, val_y


def net_force(starlist,i):
    # calls force_ij and takes the sum of all force in x & y direction
    # resulting value after for loop should be the net force on star[i]
    net_Fx, net_Fy = 0, 0
    m1, x1, y1 = starlist[i][1], starlist[i][2], starlist[i][3]
    for j in range(len(starlist)):
        if j != i:
            m2, x2, y2 = starlist[j][1], starlist[j][2], starlist[j][3]
            d = dist(x1,x2,y1,y2)
            Fx, Fy = components(force_ij(m1,m2,d),x1,x2,y1,y2)
            net_Fx += Fx
            net_Fy += Fy
    return net_Fx, net_Fy

def accel(mass,Fx,Fy):
    # calculates acceleration from net force ans mass
    # Newton's second law
    acc_x = Fx / mass
    acc_y = Fy / mass
    return acc_x, acc_y


def pos(starlist,i,ax,ay,dt):
    # calculates new position from current position, velocity, and acceleration
    _,_,x,y,vx,vy = starlist[i]
    pos_x = x + vx*dt + ax*(dt**2)/2
    pos_y = y + vy*dt + ay*(dt**2)/2
    return pos_x, pos_y


def vel(starlist,i,ax,ay,dt):
    # calculates new velocity from current velocity and acceleration
    vx,vy = starlist[i][4], starlist[i][5]
    v_x = vx + ax*dt
    v_y = vy + ay*dt
    return v_x, v_y


def starloop(starlist):
    # repeats calculation of new starlist for star[i]
    # calls net_force(), accel(), pos(), vel() within a loop for i
    new_starlist = []
    for i in range(len(starlist)):
        mass = starlist[i][1]
        Fx,Fy = net_force(starlist,i)
        ax,ay = accel(mass,Fx,Fy)
        x_new, y_new   = pos(starlist,i,ax,ay,dt)
        vx_new, vy_new = vel(starlist,i,ax,ay,dt)
        new_starlist.append([starlist[i][0],mass,x_new,y_new,vx_new,vy_new])
    return new_starlist

def time_development(initial_list,dt,t_max):
    # calls starloop() and stacks result on 'gal_hist'
    # repeats above for given time length
    gal_hist = []
    t = 0
    for k in range(int(t_max/dt)):
        if k == 0:
            starlist = starloop(initial_list)
            gal_hist.append([t,starlist])
        else:
            gal_hist.append([t,starloop(gal_hist[k-1][1])])
        t += dt
    return gal_hist
'''
def animate():
    #import visual 
    star = sphere(pos(initial condions),radius= 1,color = color.blue, make_trail = True, trail_type  = 'points', interval = 5, retain = 10 ) 
    #star velocity
    dt = 0.01 
    "  while true .... do animation" 
    # need to integrate with other functions before continuing 
    
    #edit attempt
'''

#### main ####
initial_list = data_read(csv_filename)
gal_hist = time_development(initial_list,dt,t_max)
#animate(gal_hist)
#print(*gal_hist,sep="\n")

import matplotlib.pyplot as plt
plt.figure()

for i in range(len(gal_hist)):
    finalxlist = [gal_hist[i][1][j][2] for j in range(len(gal_hist[-1][1]))]
    finalylist = [gal_hist[i][1][j][3] for j in range(len(gal_hist[-1][1]))]
    plt.scatter(finalxlist,finalylist)
    plt.grid()
    plt.show()

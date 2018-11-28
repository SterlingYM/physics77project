# constants
G = 



# function ptototypes here


def data_read(cev_filename):
    # imports csv file and returns initial condition data
    # formatted as [[star#,mass,x,y,vx,vy],[],[],...]
    import csv
    with open(cev_filename, 'rb') as csvfile:
    starlist = csv.reader(csvfile, delimiter=' ', quotechar='|')
    return  starlist


def dist(x1,x2,y1,y2):
    # calculates the distance between two sets of coordinates 
    import math
    dist = math.hypot(x2 - x1, y2 - y1)
    return dist


def force_ij(m1,m2,dist):
    # calculates the value of force between two particles given
    force_ij = (G*m1*m2)/dist
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
    def helper(j):
        if j == len(starlist):
            return net_Fx, net_Fy
        if j != i:
            m2, x2, y2 = starlist[j][1], starlist[j][2], starlist[j][3]
            d = dist(x1,x2,y1,y2)
            net_Fx, net_Fy += components(force_ij(m1,m2,d),x1,x2,y1,y2)
        helper(j+1)
    helper(0)


def acccel(mass,Fx,Fy):
    # calculates acceleration from net force ans mass
    # Newton's second law
    acc_x = Fx / mass
    acc_y = Fy / mass
    return acc_x, acc_y


def pos(starlist,i,ax,ay,dt):
    # calculates new position from current position, velocity, and acceleration
    _,_,x,y,vx,vy = starlist[i]
    pos_x = x + vx*dt + ax*(t**2)/2
    pos_y = y + vy*st + ay*(t**2)/2
    return pos_x, pos_y


def vel(starlist,i,ax,ay,dt):
    # calculates new velocity from current velocity and acceleration
    vx,vy = starlist[i][4], starlist[i][5]
    v_x = vx + ax*dt
    v_y = vy + ay*dt
    return v_x, v_y


def starloop():
    # repeats calculation of new starlist for star[i]
    # calls net_force(), accel(), pos(), vel() within a loop for i
    
    return starlist_next


def time_development():
    # calls starloop() and stacks result on 'gal_hist'
    # repeats above for given time length
    return gal_hist

def animate():
    #import visual 
    star = sphere(pos(initial condions),radius= 1,color = color.blue, make_trail = True, trail_type  = 'points', interval = 5, retain = 10 ) 
    #star velocity
    dt = 0.01 
    "  while true .... do animation" 
    # need to integrate with other functions before continuing 
    
    #edit attempt




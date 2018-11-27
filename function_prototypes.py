# function ptototypes here


def data_read():
    # import csv
    with open('title.csv', 'rb') as csvfile:
    starlist = csv.reader(csvfile, delimiter=' ', quotechar='|')
    return  starlist


def dist(x1,x2,y1,y2):
    # calculates the distance between two sets of coordinates 
    import math
    dist = math.hypot(x2 - x1, y2 - y1)
    return dist


def force_ij():
    # calculates the value of force between two particles given
    return force_ij


def components(val,x1,y1,x2,y2):
    # returns given value in component form
    # assuming the value passed is a vector on (x1,y1) to (x2,y2)
    return val_x, val_y


def net_force():
    # calls force_ij and takes the sum of all force in x & y direction
    # resulting value after for loop should be the net force on star[i]
    return net_Fx, net_Fy


def acccel():
    # calculates acceleration from net force ans mass
    # Newton's second law
    return acc_x, acc_y


def pos():
    # calculates new position from current position, velocity, and acceleration
    return pos_x, pos_y


def vel():
    # calculates new velocity from current velocity and acceleration
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




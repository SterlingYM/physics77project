import sys

year = 3600*24*365 #[sec]
kpc = 3.086 * 10**19 #[m]

def animate(saved_data):
    # animates data passed in gal_hist
    condition_data = saved_data[0]
    gal_hist       = saved_data[1]
   

    # condition data
    _,_,_,BH_mass,rho0,r_c,star_v,num_stars,num_bulge,_,_,_,_,_ = condition_data[1]

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import numpy as np
    fig = plt.figure("3D",figsize = plt.figaspect(0.7)*1.5)
    ax = fig.add_subplot(111, projection='3d')

    fig2 = plt.figure("2D",figsize = plt.figaspect(0.7)*1.5)

    for i in range(len(gal_hist)):
        ax.cla()
        # data
        xlist = [gal_hist[i][1][j][2] for j in range(len(gal_hist[-1][1]))]
        ylist = [gal_hist[i][1][j][3] for j in range(len(gal_hist[-1][1]))]
        zlist = [gal_hist[i][1][j][4] for j in range(len(gal_hist[-1][1]))]
        vxlist = [gal_hist[i][1][j][5] for j in range(len(gal_hist[-1][1]))]
        vylist = [gal_hist[i][1][j][6] for j in range(len(gal_hist[-1][1]))]
        vzlist = [gal_hist[i][1][j][7] for j in range(len(gal_hist[-1][1]))]


        # boundary
        MAX = 0.8 * 10**21
        ax.set_xlim(-MAX,MAX)
        ax.set_ylim(-MAX,MAX)
        ax.set_zlim(-MAX,MAX)
        
        # info
        info = 'initial data:\nv_avg       = {} [m/s]\nn_stars     = {} \nn_bulge     = {} \n'\
                .format(star_v,num_stars,num_bulge)+\
            'DM_rho0  = {:1.3} [kg/m^3]\nr_c            = {:1.3}[m]\nBH_mass  = {:1.3}[kg]\n'\
            .format(rho0,r_c,BH_mass)
        ax.text(MAX*0.7,MAX,MAX*1.5,info)

        # plot
        plt.figure("3D")
        ax.scatter(0,0,0,color='orange')
        ax.scatter(xlist,ylist,zlist)
        plt.title('galaxy age = {:1.2} [yr]'.format(gal_hist[i][0]/year))
        plt.pause(0.001)
    
        # plot velocity curve
        fig2.clf()
        plt.figure("2D")
        plt.xlim(0,40)
        plt.ylim(0,1500)
        plt.xlabel('radius [kpc]')
        plt.ylabel('velocity [m/s]')
        plt.title('Galaxy velocity curve')
        vel_list = np.divide(np.power(( np.power(vxlist,2) + np.power(vylist,2) + np.power(vzlist,2)),(1/2)),1000)
        r_list = np.divide(np.power( np.power(xlist,2) + np.power(ylist,2) + np.power(zlist,2),(1/2)),kpc)
        plt.scatter(r_list,vel_list)
        plt.pause(0.001)

    plt.show()


def read_data(filename):
    import pickle
    inputdata = open(filename,'rb')
    saved_data = pickle.load(inputdata)
    return saved_data


##### main #####
filename = 'gal_hist.dat'

if len(sys.argv) != 0:
    #from multiprocessing import Pool
    #saved_data = list(read_data(sys.argv[i]) for i in range(1,len(sys.argv)))
    #with Pool(len(sys.argv)) as p:
    #    p.map(animate,saved_data)
    saved_data = read_data(sys.argv[1])
    animate(saved_data)
else:
    print('Importing data from default data \'{}\'...'.format(filename))
    saved_data = read_data(filename)
    animate(saved_data)

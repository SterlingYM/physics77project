import sys

year = 3600*24*365 #[sec]

def animate(saved_data):
    # animates data passed in gal_hist
    condition_data = saved_data[0]
    gal_hist       = saved_data[1]
   

    # condition data
    _,_,_,BH_mass,rho0,r_c,star_v,num_stars,_,_,_,_ = condition_data[1]

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize = plt.figaspect(0.7)*1.5)
    ax = fig.add_subplot(111, projection='3d')

   
    for i in range(len(gal_hist)):
        ax.cla()
        # data
        xlist = [gal_hist[i][1][j][2] for j in range(len(gal_hist[-1][1]))]
        ylist = [gal_hist[i][1][j][3] for j in range(len(gal_hist[-1][1]))]
        zlist = [gal_hist[i][1][j][4] for j in range(len(gal_hist[-1][1]))]

        # boundary
        MAX = 0.8 * 10**21
        ax.set_xlim(-MAX,MAX)
        ax.set_ylim(-MAX,MAX)
        ax.set_zlim(-MAX,MAX)
        
        # info
        info = 'initial data:\nv_avg       = {} [m/s]\nn_stars     = {} \n'.format(star_v,num_stars)+\
            'DM_rho0  = {:1.3} [kg/m^3]\nr_c            = {:1.3}[m]\nBH_mass  = {:1.3}[kg]\n'\
            .format(rho0,r_c,BH_mass)
        ax.text(MAX*0.7,MAX,MAX*1.5,info)

        # plot
        ax.scatter(0,0,0,color='orange')
        ax.scatter(xlist,ylist,zlist)
        plt.title('galaxy age = {:1.2} [yr]'.format(gal_hist[i][0]/year))
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
    from multiprocessing import Pool
    saved_data = list(read_data(sys.argv[i]) for i in range(1,len(sys.argv)))
    with Pool(len(sys.argv)) as p:
        p.map(animate,saved_data)
else:
    print('Importing data from default data \'{}\'...'.format(filename))
    saved_data = read_data(filename)
    animate(saved_data)

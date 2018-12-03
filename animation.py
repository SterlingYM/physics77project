year = 3600*24*365 #[sec]

def animate(gal_hist):
    # animates data passed in gal_hist
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
        
        #plot
        ax.scatter(0,0,0,color='orange')
        ax.scatter(xlist,ylist,zlist)
        plt.title('galaxy age = {:1.2} [yr]'.format(gal_hist[i][0]/year))
        plt.pause(0.001)
        
       
    plt.show()

def read_data(filename):
    import pickle
    inputdata = open(filename,'rb')
    gal_hist = pickle.load(inputdata)
    return gal_hist

##### main #####
filename = 'gal_hist.dat'
gal_hist = read_data(filename)
animate(gal_hist)

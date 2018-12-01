import csv

num_stars = 10**2
BH_mass = 8.2 * 10**36 #[kg]
parsec = 3.086 * 10**16 #[m]
gal_radius = 25 * 10**3 * parsec #[m]
gal_disk_thickness_half = 0.15 * 10**3 * parsec #[m]
gal_bulge_radius = 0.5 * 10**3 * parsec #[m]
Msun = 10**30 #[kg]

import numpy as np

bulge_num = np.arange(num_stars)
bulge_mass = np.random.uniform(0.1*Msun,15*Msun,num_stars)
bulge_x = np.random.uniform(-1*gal_bulge_radius,gal_bulge_radius,num_stars)
bulge_y = np.random.uniform(-1*gal_bulge_radius,gal_bulge_radius,num_stars)
bulge_z = np.random.uniform(-1*gal_bulge_radius,gal_bulge_radius,num_stars)
bulge_vx = np.zeros(num_stars,dtype=float)
bulge_vy = np.zeros(num_stars,dtype=float)
bulge_vz = np.zeros(num_stars,dtype=float)


with open('Initial_Conditions_test.csv','w') as csvfile:

    csvwriter = csv.writer(csvfile, delimiter=',')
    for i in range(len(bulge_num)):
        csvwriter.writerow([bulge_num[i],bulge_mass[i],bulge_x[i],\
                bulge_y[i],bulge_z[i],bulge_vx[i],bulge_vy[i],bulge_vz[i]])

csvfile.close()

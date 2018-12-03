import csv

num_stars = 1 * 10**2
actual_num = 10 ** 9
mass_coef = actual_num / num_stars
BH_mass = 8.2 * 10**36 #[kg]
parsec = 3.086 * 10**16 #[m]
gal_radius = 25 * 10**3 * parsec #[m]
gal_disk_thickness_half = 3 * 10**3 * parsec #[m]
gal_bulge_radius = 0.5 * 10**3 * parsec #[m]
Msun = 10**30 #[kg]
star_v = 150 * 10**3 #[m/s]

import numpy as np

bulge_num = np.arange(num_stars)
bulge_mass = np.random.uniform(1*Msun*mass_coef,20*Msun*mass_coef,num_stars)
x = np.random.uniform(-1*gal_radius,gal_radius,num_stars)
y = np.random.uniform(-1*gal_radius,gal_radius,num_stars)
z = np.random.uniform(-1*gal_disk_thickness_half,gal_disk_thickness_half,num_stars)

vx = []
vy = []
vz = []
for i in range(num_stars):
    xsign = x[i]/abs(x[i])
    ysign = y[i]/abs(y[i])
    velx = np.random.uniform(star_v - 150*10**3, star_v + 150*10**3,1)
    vely = np.random.uniform(star_v - 150*10**3, star_v + 150*10**3,1)
    velz = np.random.uniform(-10*10**3,-10*10**3,1)
    vx.append(float(-1* ysign * velx))
    vy.append(float(xsign * vely))
    vz.append(float(velz))
    print(vx[i],vy[i],vz[i])
with open('Initial_Conditions_test.csv','w') as csvfile:

    csvwriter = csv.writer(csvfile, delimiter=',')
    for i in range(len(bulge_num)):
        csvwriter.writerow([bulge_num[i],bulge_mass[i],x[i],\
                y[i],z[i],vx[i],vy[i],vz[i]])

csvfile.close()

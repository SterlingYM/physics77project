import csv

num_stars = 2* 10**2
actual_num = 10 ** 10
mass_coef = actual_num / num_stars
BH_mass = 8.2 * 10**36 #[kg]
parsec = 3.086 * 10**16 #[m]
gal_radius = 25 * 10**3 * parsec #[m]
gal_disk_thickness_half = 0.15 * 10**3 * parsec #[m]
gal_bulge_radius = 0.5 * 10**3 * parsec #[m]
Msun = 1.989 * 10**30 #[kg]
star_v = 300 * 10**3 #[m/s]

import numpy as np

# position generation
bulge_num = np.arange(num_stars)
bulge_mass = np.random.uniform(1*Msun*mass_coef,20*Msun*mass_coef,num_stars)
x = np.random.triangular(-1*gal_radius,0,gal_radius,num_stars)
y = []
for i in range(len(x)):
    ylim = np.sqrt(gal_radius**2 - x[i]**2)
    y.append(*np.random.uniform(-1*ylim,ylim,1))
z = np.random.uniform(-1*gal_disk_thickness_half,gal_disk_thickness_half,num_stars)


# velocity generation
vx = []
vy = []
vz = []
for i in range(num_stars):
    d = np.sqrt(x[i]**2 + y[i]**2)
    v = np.random.uniform(star_v - 150*10**3, star_v + 250*10**3,1)
    vx.append(float(star_v * -1 * y[i] / d))
    vy.append(float(star_v * x[i] / d))
    vz.append(float(np.random.uniform(-10*10**3,-10*10**3,1)))


# data output
filename = 'Initial_Conditions_test.csv'
with open(filename,'w') as csvfile:

    csvwriter = csv.writer(csvfile, delimiter=',')
    for i in range(len(bulge_num)):
        csvwriter.writerow([bulge_num[i],bulge_mass[i],x[i],\
                y[i],z[i],vx[i],vy[i],vz[i]])
        # print(x[i],y[i],z[i])

csvfile.close()

print('Data was successfully created and saved to {}'.format(filename))

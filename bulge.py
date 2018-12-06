
bulgechar = [] #bulge characteristics 
	def bulge (gal_bulge_r,):
   radius = numpy.random.uniform(0.0,gal_bulge_r,0,(self.number_of_particles,1)) 
   theta = numpy.random.uniform(0.,1.,(self.number_of_particles,1))*pi
   phi = numpy.arccos(1-2*numpy.random.uniform(0.0,1.,(self.number_of_particles,1)))

   x = radius * numpy.sin( theta ) * numpy.cos( phi )
   y = radius * numpy.sin( theta ) * numpy.sin( phi )
   z = radius * numpy.cos( theta )
   totalmassgas = 9E18 #solar masses atomic and molecular
   m = totalmassinbulge/number_of_particles

   R  = numpy.sqrt(x**2 + y**2 + z**2)
   velc = sqrt.((G * m )/ R)


	bulgechar.append(x,y,z,m,velc)



   


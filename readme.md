note: the documentation below is written for midterm checkpoint. The structure of system has been changed for the final submission, and therefore this document is not relevant to the project anymore (Dec 2018).



# Physics 77 final project
This is a set of codes for private project.

Team member: Yukei, Tom, Savannah

Updated: 11/21/2018

Physics 77 project team 2018. All rights reserved.

--------------------

### what to do:
- N-body simulation of galaxy development(2D?)
- Animation of simulated data

### how to do:
- Take initial condition of N-particles (coordinate and velocity)
- Run gravitational simulation to calculate next position and velocity
- Loop above

- After creating data table (coordinate axis and time axis), show animation



### initial condition:
- almost uniform distribution of particles with small random fluctuation of density
- having tangential velocity against single center point (origin)


-------------------
# System Overview 
```python
# starlist: 	Table of star information at given time. 
#		Contains (1) coordinate and (2) velocity 
#		for each particle (=star).
# gal_hist: 	History of the Galaxy. 
#		Table of starlist with time axis.
```

note: function names are represented by **bold** letters

1.  **data_read()**:
	- Take in initial 'starlist'

2.  Prepare an array as gal_hist

3.  **time_development()**: for time t_next,
	1. Prepare 'starlist_next'
	2. **Starloop()**: for star[i],
		1. **net_force()**:
			- calculate force between i-j
				- **dist()**: calculate distance
				- **force_ij()**: calculate force value
			- **components()**: break force in component form
			- add all force in component form
			- get (return) net force on star[i]
		2. **accel()**:
			- get (return) acceleration of star[i] 
		3. **pos()**:
			- get (return) next coordinate from r & v & a
		4. **vel()**:
			- get (return) next velocity from v & a
		5. append star[i]'s data to 'starlist_next'
	3. vStack 'starlist_next' to 'gal_hist'
4. **animate()**:
	- show animation of position of stars over time	

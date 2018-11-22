# function ptototypes here


def data_read():
#import csv
with open('title.csv', 'rb') as csvfile:
...     starlist = csv.reader(csvfile, delimiter=' ', quotechar='|')
return  starlist

def dist(x1,x2,y1,y2):
  #import math
  #calculates the distance between two sets of coordinates 
math.hypot(x2 - x1, y2 - y1)
  
def force_ij():
 

def components():


def net_force():
  #netforce loop of ij? 

    
    

def acccel():


def pos():


def vel():


def starloop():


def time_development():


def animate():
  #import visual 
  star = sphere(pos(#initial condions),radius= 1,color = color.blue, make_trail = True, trail_type  = 'points', 
               interval = 5, retain = 10 ) 
    #star velocity
    dt = 0.01 
                
              "  while true .... do animation" 
    # need to integrate with other functions before continuing 
  




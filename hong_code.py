# Christopher Hong 
import math, sys, pygame, time
from pygame.locals import *
from datetime import datetime 

k=1
current=100 

#-------------------------------------------------------------------------------------------#

# FUNCTIONS 

# Calculating Magnetic Field
def magneticfield(x0,y0,z0,x1,y1,z1,Mx,My,Mz,Px,Py,Pz):
	r=math.sqrt((Px-Mx)**2+(Py-My)**2+(Pz-Mz)**2) 		# Distance from wire to point 
	r3inv=r**(-3)
	factor=k*current*r3inv
	dBx=factor*((y1-y0)*(Pz-Mz)-(Py-My)*(z1-z0))
	dBy=-factor*((x1-x0)*(Pz-Mz)-(Px-Mx)*(z1-z0))
	dBz=factor*((x1-x0)*(Py-My)-(Px-Mx)*(y1-y0))
	# (x0,y0,z0) is the coordinate for the initial point of the line/wire segment
	# (xl,y1,z1) is the coordinate for the final point of the line/wire segment
	# (Mx,My,Mz) is the coordinate for the midpoint of the two endpoints of the wire segment
	# (Px,Py,Pz) is the coordinate for the point in which we are finding the B-field for
	return dBx,dBy,dBz 

#-------------------------------------------------------------------------------------------#
	
# MAIN 

mult=2 

# TO COMPRESS THE SOLENOID, ADJUST num_of_points, and circle[i][1]
# i.e. num_of_points=2000 and circle[i][1]=i-1000 is more compact than 
#	   num of_points=1000 and circle[i][1]=2*(i-500)
# The x-axis bounds are from -500 to 500 
#
# thetaadd also compresses the solenoid, but it reduces the resolution of the circle
# and messes the B-field calculation up more 

# Getting the points of the circle and storing them in an array
num_of_points=10000		#Number of segments that make the polygon
num_of_points2=num_of_points/2
circle=[[0 for col in range(4)] for row in range (num_of_points+1)]

r=100 					# Adjust the radius of the circle
theta0=0
pi2=math.pi*2
thetaadd=pi2/50			# Adjust "resolution" of circle - 10 or greater recommended 
i=0
while i<=num_of_points:
	circle[i][0]=theta0
	circle[i][1]=0.5*(i-num_of_points2) 	# Adjust x-axis points
#	circle[i][1]=2*(i-num_of_points-70) 	# TO SEE FRINGING EFFECT 
	circle[1][2]=r*math.cos(theta0)			# y-axis eq for helix
	circle[i][3]=r*math.sin(theta0) 		# z-axis eq for helix 
	theta0+=thetaadd
	i+=1 

# Calculate the midpoint between two points 
midpt=[[0 for col in range(3)] for row in range (num_of_points)]
for i in range(num_of_points):
	midpt[i][0]=(circle[i][1]+circle[i+1][1])/2
	midpt[i][1]=(circle[i][2]+circle[i+1][2])/2 
	midpt[i][2]=(circle[i][3]+circle[i+1][3])/2
	
# Graphing the polygon
W,H=500,500							# Size of window
W2,H2=W/2,H/2
W20,H20=W/50,H/50 						# Adjust factor to increase/decrease # of points to calculate B-field
numofrows=357 #int(((H20)+1)*(7)) 
pygame.init() 							# Initiates pygame
screen=pygame.display.set_mode((W,H)) 	                        # Sets the window up 
pygame.display.set_caption('Magnetic Fields')                   # Caption
for i in range(num_of_points): 
	pygame.draw.aaline(screen,(75,75,75),(mult*circle[i][1]+W2,-mult*circle[i][2]+H2),(mult*circle[i+1][1]+W2,-mult*circle[i+1][2]+H2),2)

pygame.display.flip()
# Setting an array with all the points to be tested 
P=[[0 for col in range(3)] for row in range(numofrows)]         # numofrows SHOULD BE 357!!!

i=0
for Pi in range(int(-2*W/25-13),int(2*W/25),int(H20+5)):
	for Pj in range(int(-H2),int(H2+1),int(H20)):
		print("starting i=", i)
		P[i][0]=Pi
		print("Pi placed successfully: ", P[i][0])
		P[i][1]=Pj
		print("Pj placed successfully: ", P[i][1])
		#print(P)
		i+=1


# Calculating the B Field at point P[i]
dBtotal=[[0 for col in range(3)] for row in range(numofrows)] 

# initiate loop here
for Pi in range(numofrows):
	for i in range(num_of_points):
		i1=i+1
		x0=circle[i][1]
		y0=circle[i][2]
		z0=circle[i][3]
		x1=circle[i1][1]
		y1=circle[i1][2]
		z1=circle[i1][3]
		Mx=midpt[i][0]
		My=midpt[i][1]
		Mz=midpt[i][2]
		Px=P[Pi][0]
		Py=P[Pi][1]
		Pz=P[Pi][2]
		dBx,dBy,dBz=magneticfield(x0,y0,z0,x1,y1,z1,Mx,My,Mz,Px,Py,Pz) # B-field function 
		dBtotal[Pi][0]+=dBx 										   # Calculate total x-comp. B-field 
		dBtotal[Pi][1]+=dBy 										   # Calculate total y-comp. B-field 
		dBtotal[Pi][2]+=dBz 										   # Calculate total z-comp. B-field 

		
# Graphing the B-field
magB=[[0 for col in range(1)] for row in range(numofrows)]
magBmax=10**(-10)
magBmin=10**10
for Pi in range(numofrows):
	magB[Pi]=math.sqrt(dBtotal[Pi][0]**2+dBtotal[Pi][1]**2+dBtotal[Pi][2]**2) # magnitude of B field at point
	if magB[Pi] < magBmin:		# get max B magnitude
		magBmin = magB[Pi]
	if magB[Pi] > magBmax: 		# get min B magnitude
		magBmax = magB[Pi]
		

# Set scale for intensity of color
rangeofmagB=magBmax-magBmin
colorscale=190/rangeofmagB

# Slope of B-field at point P
slopeBfield=[[0 for col in range(1)] for row in range(numofrows)]
 
for Pi in range(numofrows):
	if dBtotal[Pi][0] == 0:
		slopeBfield[Pi]=10
	if dBtotal[Pi][0] != 0:
		slopeBfield[Pi]=dBtotal[Pi][1]/dBtotal[Pi][0] 


# Graph B-field
for Pi in range(numofrows): 
	color=magB[Pi]*colorscale
	pygame.draw.aaline(screen,(0,color+55,0),(mult*(P[Pi][0]-2.5*dBtotal[Pi][0]/magB[Pi])+W2,-(mult*(P[Pi][1]-2.5*dBtotal[Pi][1]/magB[Pi]))+H2),(mult*(P[Pi][0]+2.5*dBtotal[Pi][0]/magB[Pi])+W2,-(mult*(P[Pi][1]+2.5*dBtotal[Pi][1]/magB[Pi]))+H2),10)
	pygame.draw.aaline(screen,(0,color+55,0),(mult*(P[Pi][0]-2.5*dBtotal[Pi][0]/magB[Pi])+W2,-(mult*(P[Pi][1]-2.5*dBtotal[Pi][1]/magB[Pi]))+H2-1),(mult*(P[Pi][0]+2.5*dBtotal[Pi][0]/magB[Pi])+W2,-(mult*(P[Pi][1]+2.5*dBtotal[Pi][1]/magB[Pi]))+H2-1),10)
	pygame.display.flip()
	
running=True 

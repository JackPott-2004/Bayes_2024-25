#A rocket path prediction simulation
#the necessary imports
import matplotlib.pyplot as plt
import numpy as np
import json

timestep = 0.052 # Isaac said just over 50ms

orientation_x = []
orientation_y = []
orientation_z = []
linAccel_x = []
linAccel_y = []
linAccel_z = []

#read the data from the CSV file
with open('combined.json') as file:
    for line in file:
        
        line = line.replace('"altitude": }', '"altitude": null}')
        rocket_parameters = json.loads(line) #json.load(file)
        orientation_x.append(rocket_parameters["bno055"]["orientation"]["x"])
        orientation_y.append(rocket_parameters['bno055']["orientation"]["y"])
        orientation_z.append(rocket_parameters['bno055']["orientation"]["z"])
        linAccel_x.append(rocket_parameters['bno055']["linear_acceleration"]["x"])
        linAccel_y.append(rocket_parameters['bno055']["linear_acceleration"]["y"])
        linAccel_z.append(rocket_parameters['bno055']["linear_acceleration"]["z"])

size = len(orientation_x)

velo_x = []
velo_y = []
velo_z = []

posi_x = []
posi_y = []
posi_z = []

a_x = []
a_y = []
a_z = []

def changeInSpeed(acceleration):
    v = acceleration * timestep 
    return v

def changeInPosition(velocity, acceleration):
    p = velocity * timestep + 0.5 * acceleration * ((timestep)**2)
    return p                                                                                                                                                                     

def velocitySum(ax,ay,az):
    speed_x = 0.0
    speed_y = 0.0
    speed_z = 0.0
    for i in range(size):

        speed_x += changeInSpeed(ax[i])
        speed_y += changeInSpeed(ay[i])
        speed_z += changeInSpeed(az[i])
        velo_x.append(speed_x)
        velo_y.append(speed_y)
        velo_z.append(speed_z)
        
    return velo_x, velo_y, velo_z

def positionSum(vx,vy,vz,ax,ay,az):
    r_x = 0.0
    r_y = 0.0
    r_z = 0.0
    for i in range(size):
        r_x += changeInPosition(vx[i],ax[i])
        r_y += changeInPosition(vy[i],ay[i])
        r_z += changeInPosition(vz[i],az[i])
        posi_x.append(r_x)
        posi_y.append(r_y)
        posi_z.append(r_z)
        
    return posi_x, posi_y, posi_z

def plotter_3D(x,y,z):
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot3D(x,y,z, linewidth=5)
    """
    ax.set_xlim(0, 10000)  # Adjust X-axis range
    ax.set_ylim(0, 10000)  # Adjust Y-axis rang
    ax.set_zlim(0, 10000)  # Adjust Z-axis range
    """
    ax.set_xlabel('Z-axis') 
    ax.set_ylabel('Y-axis')  
    ax.set_zlabel('X-axis')  
    plt.show()

def main():
    ax, ay, az = -1 * np.array(linAccel_z), -1 *  np.array(linAccel_y), -1*np.array(linAccel_x)
    vx, vy, vz    = velocitySum(ax, ay, az)
    px, py, pz    = positionSum(vx, vy, vz, ax, ay, az)
    plotter_3D(px,py,pz)
        

main()
"""
time = numbers = list(range(1, len(orientation_x) + 1))
plt.plot(time, orientation_y )
plt.show()"

"""
#A rocket path prediction simulation
"""This program aims to predict the path of the rocket which the bayes "PCB" will be flying on to gather test data
It does this by figuring out the acceleration in the grounds frame by transforming the acceleration in the rockets frame through rotation matrices,
as we have data on the rockets
"""
#the necessary imports
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation as R
import json

timestep = 0.1

orientation_x = []
orientation_y = []
orientation_z = []
linAccel_x = []
linAccel_y = []
linAccel_z = []
accelerometer_x = []
accelerometer_y = []
accelerometer_z = []
linAcc = [linAccel_x, linAccel_y, linAccel_z]


#read the data from the CSV file
with open('combined.json') as file:
    for line in file:
        line = line.strip()

        if not line:
            continue

        try:
            line = line.replace('"altitude": }', '"altitude": null}')
            rocket_parameters = json.loads(line) #json.load(file)
            orientation_x.append(rocket_parameters["bno055"]["orientation"]["x"])
            orientation_y.append(rocket_parameters['bno055']["orientation"]["y"])
            orientation_z.append(rocket_parameters['bno055']["orientation"]["z"])
            linAccel_x.append(rocket_parameters['bno055']["linear_acceleration"]["x"])
            linAccel_y.append(rocket_parameters['bno055']["linear_acceleration"]["y"])
            linAccel_z.append(rocket_parameters['bno055']["linear_acceleration"]["z"])
            accelerometer_x.append(rocket_parameters['bno055']["accelerometer"]["x"])
            accelerometer_y.append(rocket_parameters['bno055']["accelerometer"]["y"])
            accelerometer_z.append(rocket_parameters['bno055']["accelerometer"]["z"])

        except json.JSONDecodeError as e:
            print(f"JSON Decode Error: {e} in line: {line}")  # Debugging

velo_x = []
velo_y = []
velo_z = []

posi_x = []
posi_y = []
posi_z = []

a_x = []
a_y = []
a_z = []


def placeholder(linAccel_x,linAccel_y, linAccel_z):
    for i in range(len(linAccel_x)):
        theta = float(orientation_x[i])
        phi   = float(orientation_y[i])
        
        accel_x = float(linAccel_x[i])
        accel_y = float(linAccel_y[i])
        accel_z = float(linAccel_z[i])

        rotation_pitch = np.array([[1, 0, 0] ,
                                   [0, np.cos(theta), np.sin(theta)],
                                   [0, -1 * np.sin(theta), np.cos(theta)]])
        
        rotation_yaw   = np.array([[np.cos(phi), 0, np.sin(phi)],
                                   [0, 1, 0],
                                   [-1 * np.sin(phi), 0 , np.cos(phi)]])

        rotation_total = np.dot(rotation_pitch, rotation_yaw)
        
        a_prime = np.array([[accel_x], [accel_y], [accel_z]])
        a = np.dot(rotation_total, a_prime)
        a_x.append(float(a[0][0]))
        a_y.append(float(a[1][0]))
        a_z.append(float(a[2][0]))
    return a_x, a_y, a_z

def changeInSpeed(acceleration):
    """delta_v = at"""
    v = acceleration * timestep 
    return v

def changeInPosition(velocity, acceleration):
    """Uses the s = ut + 0.5at^2"""
    p = velocity * timestep + 0.5 * acceleration * ((timestep)**2)
    return p                                                                                                                                                                     


def velocitySum(ax,ay,az):
    """This creates an array for the speed in each grou axis direction"""
    speed_x = 0.0
    speed_y = 0.0
    speed_z = 0.0
    for i in range(len(ax)):
        #i swapped these blocks around
        velo_x.append(speed_x)
        velo_y.append(speed_y)
        velo_z.append(speed_z)
        speed_x += changeInSpeed(ax[i])
        speed_y += changeInSpeed(ay[i])
        speed_z += changeInSpeed(az[i])
       
    return velo_x, velo_y, velo_z

def positionSum(vx,vy,vz,ax,ay,az):
    r_x = 0.0
    r_y = 0.0
    r_z = 0.0
    for i in range(len(ax)):
        #i swapped these blocks around
        posi_x.append(r_x)
        posi_y.append(r_y)
        posi_z.append(r_z)
        r_x += changeInPosition(vx[i],ax[i])
        r_y += changeInPosition(vy[i],ay[i])
        r_z += changeInPosition(vz[i],az[i])
        
    return posi_x, posi_y, posi_z


def plotter_3D(x,y,z):
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111, projection='3d')
    
    ax.plot3D(x,y,z, linewidth=5)
    ax.set_xlabel('X-axis') 
    ax.set_ylabel('Y-axis')  
    ax.set_zlabel('Z-axis')  
    plt.show()

def main():
    ax, ay, az = placeholder(linAccel_x, linAccel_y, linAccel_z)
    vx, vy, vz    = velocitySum(linAccel_x, linAccel_y, linAccel_z)
    px, py, pz    = positionSum(vx, vy, vz, ax, ay, az)
    plotter_3D(px,py,pz)
        
main()


"""

"""




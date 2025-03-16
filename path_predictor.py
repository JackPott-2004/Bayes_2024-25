#A rocket path prediction simulation
#the necessary imports
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation as R
import json

timestep = 0.02

orientation_x = []
orientation_y = []
orientation_z = []
linAccel_x = []
linAccel_y = []
linAccel_z = []

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

        except json.JSONDecodeError as e:
            print(f"JSON Decode Error: {e} in line: {line}")  # Debugging

    #return orientation_x, orientation_y, orientation_z, linAccel_x, linAccel_y, linAccel_z


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


def placeholder(ox,oy,linAccel_x,linAccel_y, linAccel_z,method):
    #These are no longer used as the BNO does it for me, what a waste of time 
    for i in range(size):
        if method == 1:
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
            
            return 
         
        if method == 2:
        
            accel_array = np.array([linAccel_x[i], linAccel_y[i], linAccel_z[i]])

            theta = ox[i]
            phi   = oy[i]

            axis_x = np.array([1,0,0])
            axis_y = np.array([0,1,0])

            q_x = R.from_rotvec(theta * axis_x)
            q_y = R.from_rotvec(phi * axis_y)

            total_rotation = q_x * q_y
            a = total_rotation.apply(accel_array) 
                
            #a_prime = np.array([linAccel_x, linAccel_y, linAccel_z])
            a =  total_rotation.apply(accel_array)
            a_x.append(float(a[0]))
            a_y.append(float(a[1]))
            a_z.append(float(a[2]))
        
    return a_x, a_y, a_z

def changeInSpeed(acceleration):
    """delta_v = at"""
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
    ax.set_xlabel('X-axis') 
    ax.set_ylabel('Y-axis')  
    ax.set_zlabel('Z-axis')  
    plt.show()

def main():
    method = 2
    #ax, ay, az = placeholder(orientation_x, orientation_y, -1* np.array(linAccel_x), -1 *np.array(linAccel_y), np.array(linAccel_z), method)
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
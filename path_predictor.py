#A rocket path prediction simulation
"""This program aims to predict the path of the rocket which the bayes "PCB" will be flying on to gather test data
It does this by figuring out the acceleration in the grounds frame by transforming the acceleration in the rockets frame through rotation matrices,
as we have data on the rockets
"""
#the necessary imports
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np

file  = "f.txt"

timestep = 1

columns = [[] for _ in range(29)]

#read the data from the CSV file
with open(file,"r") as file:
    lines = file.readlines()
    line_count = len(lines)
    for line in lines:
        values = line.strip().split(',')
        for i in range(len(values)):
            columns[i].append(values[i])

#BNO055
times = columns[0]
orientation_x = columns[1]
orientation_y = columns[2]
orientation_z = columns[3]
angVelocity_x = columns[4]
angVelocity_y = columns[5]
angVelocity_z = columns[6]
linAccel_x    = columns[7]
linAccel_y    = columns[8]
linAccel_z    = columns[9]
Magnetometer_x = columns[10]
Magnetometer_y = columns[11]
Magnetometer_z = columns[12]
accelerometer_x = columns[13]
accelerometer_y = columns[14]
accelerometer_z = columns[15]
gravVector_x = columns[16]
gravVector_y = columns[17]
gravVector_z = columns[18]
calibration_system = columns[19]
calibration_gyro = columns[20]
calibration_accel = columns[21]
calibration_mag = columns[22]
boardTemp = columns[23]
#
sensorData_temp = columns[24]
sensorData_pres = columns[25]
#BMP280
BMP280Temp = columns[26]
BMP280Pres = columns[27]
BMP280Alt  = columns[28]

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


def placeholder():
    for i in range(size):
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
    for i in range(size):
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
    for i in range(size):
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

    ax.set_xlim(0, 50)  # Adjust X-axis range
    ax.set_ylim(0, 50)  # Adjust Y-axis rang
    ax.set_zlim(0, 125)  # Adjust Z-axis range
    ax.set_xlabel('X-axis') 
    ax.set_ylabel('Y-axis')  
    ax.set_zlabel('Z-axis')  
    plt.show()

def main():
    ax, ay, az = placeholder()
    vx, vy, vz    = velocitySum(ax, ay, az)
    px, py, pz    = positionSum(vx, vy, vz, ax, ay, az)
    plotter_3D(px,py,pz)
        
main()


"""
NOTES:
the structure of this code is rather poor IMO, could be worked upon in the future and probably easily made OOP'esque
"""
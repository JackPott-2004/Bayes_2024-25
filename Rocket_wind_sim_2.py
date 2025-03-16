"""
Plan:
use the one class for one thing idea
make class for the rocket
make class for the canards
make a class to compute the simulation
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation as R
import math


class Rocket():
    
    def __init__(self):
        self.THRUST = 300.0
        self.LENGTH = 1.5
        self.MASS = 10
        self.GRAVITY = -9.81
        self.TIMESTEP = 0.1
        self.DISTANCE = -0.2
        self.INERTIA =  (1/12) * self.LENGTH**2 * self.MASS*0.05
        self.DRAG_CONSTANT = 0.3
        self.RADIUS = 0.1
        self.AREA = self.RADIUS * 2 * self.LENGTH
        

    def windTorqueCalc(self,orientations,RWS): 
        """This method calculates the torque on the rocket"""
    
        #As the rockets orientation changes there will be less surface area in the direction of the wind to be hit
        crossSectionalArea_x = self.AREA * (np.cos(orientations[0])) 
        crossSectionalArea_y = self.AREA * (np.cos(orientations[1]))
        #crossSectionalArea_z = self.AREA * (np.cos(self.ORIENTATIONS[2]))
        
        force_x =  ((RWS[0])**2) * crossSectionalArea_x * self.DRAG_CONSTANT * np.sign(RWS[0])
        force_y =  ((RWS[1])**2) * crossSectionalArea_y * self.DRAG_CONSTANT * np.sign(RWS[1])
        #force_z =  ((self.RWS[2])**2) * crossSectionalArea_z * self.DRAG_CONSTANT * np.sign(self.RWS[2])
        force_vector = np.array([force_x, force_y,0])
        moment_arm =np.array([0, 0, -self.DISTANCE])
        torque = np.cross(moment_arm, force_vector)
        torque_x =  torque[0]
        torque_y =  torque[1]
        return torque_x, torque_y
    
class Canard():
    
    def __init__(self):
        self.CANARD_AREA = 0.4
        self.CANARD_DISTANCE = 0.2 
        self.CANARD_RADIUS = 0.2
        self.CANARD_FORCE = [0,0]
        self.ORIENTATIONS = [0,0] 
    
    
    def canardTorqueCalc(self,RWS,orientations):
        """ Currently experiemnting with these RWS being untransformed"""
        area = self.CANARD_AREA
        area_X = area * np.abs(np.cos(orientations[1] - self.ORIENTATIONS[0]))
        area_Y = area * np.abs(np.cos(orientations[1] - self.ORIENTATIONS[1]))
        area_ZX = area * np.abs(np.sin(orientations[1] - self.ORIENTATIONS[0]))
        area_ZY = area * np.abs(np.sin(orientations[1] - self.ORIENTATIONS[1]))

        forceOnXCanard_LAT = area_X * RWS[0]**2 * np.sign(RWS[0])
        forceOnYCanard_LAT = area_Y * RWS[1]**2 * np.sign(RWS[1])
        forceOnXCanard_VERT = area_ZX * RWS[2]**2 * np.sign(RWS[2])
        forceOnYCanard_VERT = area_ZY * RWS[2]**2 * np.sign(RWS[2])
        
        force_vector = np.array([forceOnXCanard_LAT, forceOnYCanard_LAT, (forceOnYCanard_VERT + forceOnXCanard_VERT)])
        moment_arm =np.array([0, 0, self.CANARD_DISTANCE])
        torque = np.cross(moment_arm, force_vector)
        torque_x = torque[0]
        torque_y = torque[1]
        return torque_x, -torque_y
    

    def proportional(self,error):
        GAIN = 0
        per_x = error[0]
        per_y = error[1]
        self.ORIENTATIONS[0] = -(np.pi / 4) * per_x * GAIN
        self.ORIENTATIONS[1] = -(np.pi / 4) * per_y * GAIN
        return self
    
class execute():

    def __init__(self):
        self.POSITIONS = [0,0,0]
        self.VELOCITIES = [0,0,0]
        self.ORIENTATIONS = [0,0]
        self.ANGULAR_VELOCITIES = [0,0]
        self.TIMESTEP = 0.05
        self.rocket = Rocket()
        self.canard = Canard()
    
    def inputs(self):
        WS_X = 50.0 #input("What is the wind speed coming in the x direction: ")
        WS_Y = 0 #input("What is the wind speed coming in the y direction: ")
        WS_Z = 0.0
        self.DURATION = 30 #input("How long should the rocket fly for: ")
        self.WS = [WS_X, WS_Y, WS_Z]

    def relativeWindSpeedCalc(self):
        RWS_X = -self.VELOCITIES[0] + self.WS[0]
        RWS_Y = -self.VELOCITIES[1] + self.WS[1]
        RWS_Z = -self.VELOCITIES[2] + self.WS[2]
        self.RWS = [RWS_X, RWS_Y, RWS_Z]
        #print(f"The RWS in x,y,z are {self.RWS}")
        return self

    def rotateToRocketsAxis(self, array):
        theta = self.ORIENTATIONS[0]
        phi = self.ORIENTATIONS[1]
        if theta == 0 or phi == 0:
            return array

        EPSILON = 1e-6  # Small value to prevent zero rotation

        # If angles are too small, use the identity rotation
        if abs(theta) < EPSILON and abs(phi) < EPSILON:
            return array  # No rotation needed

        axis_x = np.array([1,0,0])
        axis_y = np.array([0,1,0])

        # Ensure no zero-norm quaternion
        q_x = R.from_rotvec((theta if abs(theta) >= EPSILON else EPSILON) * axis_x)
        q_y = R.from_rotvec((phi if abs(phi) >= EPSILON else EPSILON) * axis_y)

        total_rotation = q_y * q_x
        return total_rotation.apply(array)
    

    def windForce(self):
        theta = self.ORIENTATIONS[0]
        phi = self.ORIENTATIONS[1]

        EPSILON = 1e-6
        if abs(theta) < EPSILON:
            theta = EPSILON
        if abs(phi) < EPSILON:
            phi = EPSILON

        crossSectionalArea_X = self.rocket.LENGTH * self.rocket.RADIUS  * (np.cos(theta)) 
        crossSectionalArea_Y = self.rocket.LENGTH * self.rocket.RADIUS * (np.cos(phi))
        crossSectionalArea_Z = self.rocket.RADIUS**2 * np.pi * 0.5
        
        
        self.relativeWindSpeedCalc()
        force_X =  -((self.RWS[0])**2) * crossSectionalArea_X * self.rocket.DRAG_CONSTANT
        force_Y =  -((self.RWS[1])**2) * crossSectionalArea_Y * self.rocket.DRAG_CONSTANT  
        force_Z =  -((self.RWS[2])**2) * crossSectionalArea_Z * self.rocket.DRAG_CONSTANT  #includes nose cone drag coeffiecient
        
        wind_force = [force_X, force_Y, force_Z]

        return wind_force

    def rotatedThrust(self):
        accel_array = [0,0,self.rocket.THRUST]
        rotated = self.rotateToRocketsAxis(accel_array)
        return rotated

    def changeInVelocity(self):
        """This need to be redone"""
        wind_force = self.windForce()
        rotated_thrust = self.rotatedThrust()
        #print(wind_force)
        self.VELOCITIES[0] += (rotated_thrust[0] + -wind_force[0]) / self.rocket.MASS * self.TIMESTEP
        self.VELOCITIES[1] += (rotated_thrust[1] + wind_force[1]) / self.rocket.MASS * self.TIMESTEP
        self.VELOCITIES[2] += (rotated_thrust[2] + wind_force[2] + self.rocket.GRAVITY * self.rocket.MASS) / self.rocket.MASS * self.TIMESTEP
        return self
    
    def changeInPosition(self):
        self.POSITIONS[0] += self.VELOCITIES[0] * self.TIMESTEP
        self.POSITIONS[1] += self.VELOCITIES[1] * self.TIMESTEP
        self.POSITIONS[2] += self.VELOCITIES[2] * self.TIMESTEP
        return self
    
    def changeInAngularVelocity(self, torque_x, torque_y, canard_torque_x, canard_torque_y,):
        """This method finds the change (in radians per second) of the rockets angular velocities"""
        #print(f"the y torques are:{(torque_y + canard_torque_y)}")
        #print(f"the x torques are:{(torque_x + canard_torque_x)}")
        self.ANGULAR_VELOCITIES[0] = self.TIMESTEP * (torque_y + canard_torque_y) / self.rocket.INERTIA # here too
        self.ANGULAR_VELOCITIES[1] = self.TIMESTEP * (torque_x + canard_torque_x)/  self.rocket.INERTIA
        return self
        

    def changeInOrientation(self):
        """This method finds the change (in radians) of the rockets orientation"""
        self.ORIENTATIONS[0] += self.ANGULAR_VELOCITIES[0] * self.TIMESTEP
        self.ORIENTATIONS[1] += self.ANGULAR_VELOCITIES[1] * self.TIMESTEP
        return self

    def error(self):
        actual = [self.ORIENTATIONS[0], self.ORIENTATIONS[1]]
        normalised_x = -actual[0]/ (2 * np.pi)
        normalised_y = -actual[1]/ (2 * np.pi)
        return [normalised_x, normalised_y ]

    def simulate(self):
        n = 0
        Ps = [[],[],[]]
        Os = [[],[]]
        Vs = [[],[],[]]
        AVs = [[],[]]
        time = 0
        while time < self.DURATION:
            Ps[0].append(self.POSITIONS[0])
            Ps[1].append(self.POSITIONS[1])
            Ps[2].append(self.POSITIONS[2])
            Os[0].append(self.ORIENTATIONS[0])
            print(self.ORIENTATIONS[0])
            Os[1].append(self.ORIENTATIONS[1])       
            Vs[0].append(self.VELOCITIES[0]) 
            Vs[1].append(self.VELOCITIES[1]) 
            Vs[2].append(self.VELOCITIES[2]) 
            AVs[0].append(self.ANGULAR_VELOCITIES[0]) 
            AVs[1].append(self.ANGULAR_VELOCITIES[1]) 
            

            self.changeInPosition()
            self.changeInOrientation()
            self.changeInVelocity()
            #print(self.VELOCITIES)
            #print(self.ORIENTATIONS)
             
            #self.relativeWindSpeedCalc()
            RWS = [self.RWS[0], self.RWS[1], self.RWS[2]]
            RWS_transformed = self.rotateToRocketsAxis(RWS)
            RTx, RTy = self.rocket.windTorqueCalc(self.ORIENTATIONS, RWS)
            CTx, CTy = self.canard.canardTorqueCalc(RWS_transformed, self.ORIENTATIONS)
            self.changeInAngularVelocity(RTx, RTy,CTx, CTy)
            error = self.error()
            #self.canard.proportional(error)

            time += self.TIMESTEP
            if time > 10: 
                self.rocket.THRUST = 0 
               
        return Ps, Os, Vs, AVs 

"""FUNCTIONS """
def plotter_3D(x,y,z):
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot3D(x,y,z, marker='o')
    ax.set_xlabel('X-axis') 
    ax.set_ylabel('Y-axis')  
    ax.set_zlabel('Z-axis')  
    """ax.set_xlim3d(0,1000)
    ax.set_ylim3d(0,1000)
    ax.set_zlim3d(0,1000)"""
    plt.show()

def plot_2d(x,y):
    plt.plot(x,y, color='b', label='Line 1')
    plt.show()

def main():
    system = execute()
    system.inputs()
    Ps, Os, Vs, AVs = system.simulate()
    length = len(Ps[0])
    times = np.linspace(0,10,int(system.DURATION/system.TIMESTEP))
    plotter_3D(Ps[0], Ps[1], Ps[2])
    #plotter_3D(times, Os[0], Os[1])
    #plot_2d(times,Os[0])

main()
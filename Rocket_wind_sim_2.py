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


class Rocket():
    
    def __init__(self):
        self.THRUST = 30.0
        self.LENGTH = 1.0
        self.MASS = 5
        self.GRAVITY = -9.81
        self.AREA = 0.2
        self.TIMESTEP = 0.1
        self.DISTANCE = 0.2
        self.INERTIA = 0.1
        self.DRAG_CONSTANT = 0.3
        self.RADIUS = 0.1

    def windTorqueCalc(self): 
        """This method calculates the torque on the rocket"""
    
        #As the rockets orientation changes there will be less surface area in the direction of the wind to be hit
        crossSectionalArea_x = self.AREA * (np.cos(self.ORIENTATIONS[0])) 
        crossSectionalArea_y = self.AREA * (np.cos(self.ORIENTATIONS[1]))
        #crossSectionalArea_z = self.AREA * (np.cos(self.ORIENTATIONS[2]))
        
        self.relativeWindSpeedCalc()
        force_x =  ((self.RWS[0])**2) * crossSectionalArea_x * self.DRAG_CONSTANT * np.sign(self.RWS[0])
        force_y =  ((self.RWS[1])**2) * crossSectionalArea_y * self.DRAG_CONSTANT * np.sign(self.RWS[1])
        #force_z =  ((self.RWS[2])**2) * crossSectionalArea_z * self.DRAG_CONSTANT * np.sign(self.RWS[2])
        force_vector = np.array([force_x, force_y,0])
        moment_arm =np.array([0, 0, -self.DISTANCE])
        torque = np.cross(moment_arm, force_vector)
        print(torque)
        torque_x =  torque[0]
        torque_y =  torque[1]

        return torque_x, torque_y
    
    
    


class Canard():
    
    def __init__(self):
        self.CANARD_AREA = 1
        self.CANARD_DISTANCE = 0.2 
        self.CANARD_RADIUS = 0.2
        self.CANARD_VECTORS = [[self.CANARD_RADIUS,0,self.CANARD_DISTANCE], # right x one
                               [0,self.CANARD_RADIUS,self.CANARD_DISTANCE], # the top y one 
                               [-self.CANARD_RADIUS,0,self.CANARD_DISTANCE], # the left x one 
                               [0,self.CANARD_RADIUS,self.CANARD_DISTANCE]] # the bottom y one
        self.CANARD_FORCE = [0,0]
        self.ORIENTATIONS = [0,0] 

    

    def torqueCalc():
        p=p


class execute():

    def __init__(self):
        self.POSITIONS = [0,0,0]
        self.VELOCITIES = [0,0,0]
        self.ORIENTATIONS = [0,0,0]
        self.TIMESTEP = 0.05
        rocket = Rocket()
        canard = Canard()
    
    def inputs(self):
        WS_X = 10 #input("What is the wind speed coming in the x direction: ")
        WS_Y = 0 #input("What is the wind speed coming in the y direction: ")
        WS_Z = 0
        self.DURATION = 10 #input("How long should the rocket fly for: ")
        self.WS = [WS_X, WS_Y, WS_Z]

    def relativeWindSpeedCalc(self):
        RWS_X = self.velocity[0] - self.WS[0]
        RWS_Y = self.velocity[1] - self.WS[1]
        RWS_Z = self.velocity[2] - self.WS[2]
        self.RWS = [RWS_X, RWS_Y, RWS_Z]

    def windForce(self,rocket):
        accel_array = [0,0,rocket.THRUST]

        theta = self.ORIENTATIONS[0]
        phi   = self.ORIENTATIONS[1]

        axis_x = np.array([1,0,0])
        axis_y = np.array([0,1,0])

        q_x = R.from_rotvec(theta * axis_x)
        q_y = R.from_rotvec(phi * axis_y)

        total_rotation = q_x * q_y
        a = total_rotation.apply(accel_array)

        crossSectionalArea_X = rocket.LENGTH * rocket.RADIUS  * (np.cos(theta)) 
        crossSectionalArea_Y = rocket.LENGTH * rocket.RADIUS * (np.cos(phi))
        crossSectionalArea_Z = rocket.RADIUS**2 * np.pi() * 0.5
        
        
        self.relativeWindSpeedCalc()
        force_X =  ((self.RWS[0])**2) * crossSectionalArea_X * rocket.DRAG_CONSTANT
        force_Y =  ((self.RWS[1])**2) * crossSectionalArea_Y * rocket.DRAG_CONSTANT
        force_Z =  ((self.RWS[2])**2) * crossSectionalArea_Z * rocket.DRAG_CONSTANT * 0.5 #includes nose cone drag coeffiecient
        
        force_array = np.array([force_X, force_Y, force_Z])
        wind_force_rotated = total_rotation.apply(force_array)

        return a,wind_force_rotated
    


    def totalTorque(self,rocket,canard, wind_force_rotated):
        p


    def changeInVelocity(self,rocket):
        a, wind_force_rotated = self.windForce(rocket)

        self.VELOCITIES[0] += (a[0] + wind_force_rotated[0] / self.MASS) * self.TIMESTEP
        self.VELOCITIES[1] += (a[1] + wind_force_rotated[1] / self.MASS) * self.TIMESTEP
        self.VELOCITIES[2] += (a[2] + wind_force_rotated[2] + rocket.GRAVITY) * self.TIMESTEP
        return self
    
    def changeInPosition(self):
        self.POSITIONS[0] += self.VELOCITIES[0] * self.TIMESTEP
        self.POSITIONS[1] += self.VELOCITIES[1] * self.TIMESTEP
        self.POSITIONS[2] += self.VELOCITIES[2] * self.TIMESTEP
        return self
    
    def changeInAngularVelocity(self, torque_x, torque_y):
        """This method finds the change (in radians per second) of the rockets angular velocities"""
        self.ANGULAR_VELOCITY[0] = self.TIMESTEP * (torque_x) / self.INERTIA # here too
        self.ANGULAR_VELOCITY[1] = self.TIMESTEP * (torque_y)/ self.INERTIA
        return self

    def changeInOrientation(self):
        """This method finds the change (in radians) of the rockets oreintation"""
        self.ORIENTATIONS[0] += self.ANGULAR_VELOCITY_X * self.TIMESTEP
        self.ORIENTATIONS[1] += self.ANGULAR_VELOCITY_Y * self.TIMESTEP
        return self

    def error(self):
        p=p

    def simulate(self):
        Ps = [[],[],[]]
        Os = [[],[],[]]
        time = 0
        while time < self.DURATION:
            print("here")
            Ps[0].append(self.POSITIONS[0])
            Ps[1].append(self.POSITIONS[1])
            Ps[2].append(self.POSITIONS[2])
            Os[0].append(self.ORIENTATIONS[0])
            Os[1].append(self.ORIENTATIONS[1]) 
            Os[2].append(self.ORIENTATIONS[2]) 

            self.changeInPosition()
            self.changeInOrientation()
            self.changeInVelocity()
            self.changeInAngularVelocity()

            time += self.TIMESTEP
        
        return Ps, Os

"""FUNCTIONS """
def plotter_3D(x,y,z):
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot3D(x,y,z, linewidth=5)
    ax.set_xlabel('X-axis') 
    ax.set_ylabel('Y-axis')  
    ax.set_zlabel('Z-axis')  
    plt.show()

def plot_2d(x,y):
    plt.plot(x,y, marker='o', linestyle='-', color='b', label='Line 1')
    plt.show()

def main():
    system = execute()
    system.inputs()
    Ps, Os = system.simulate()
    plotter_3D(Ps[0], Ps[1], Ps[2])

main()
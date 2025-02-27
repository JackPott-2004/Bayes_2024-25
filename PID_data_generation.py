#Rocket Sim with wind

"""
The rocket should never be more than ___ degrees from the vertical so

FINAL RESULT
you will be prompted how long you'd like the rocket to fly for
you will be prompted for the wind speed in x direction
you will be prompted for the wind speed in y direction
the rockets trajectory will be plotted, it will show it being cocked into the wind
optional - put a delay in when wind starts, variable wind speed with height


"""
 
#The imports
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np



class Rocket():
    
    def __init__(self):
        self.LENGTH = 1.0
        self.MASS = 5
        self.AREA = 0.2
        self.TIMESTEP = 0.01
        #From centre of pressure to the centre of gravity,Should cock into positive wind values
        self.DISTANCE = 0.2 
        #initially vertical
        self.ORIENTATION_X = 0.0 #radians
        self.ORIENTATION_Y = 0.0 #radians
        self.ORIENTATION_Z = 0.0 #radians
        self.ANGULAR_VELOCITY_X = 0.0
        self.ANGULAR_VELOCITY_Y = 0.0
        self.ANGULAR_VELOCITY_Z = 0.0
        self.VELOCITY_X = 0.0
        self.VELOCITY_Y = 0.0
        self.VELOCITY_Z = 0.0
        self.POSITION_X = 0.0
        self.POSITION_Y = 0.0
        self.POSITION_Z = 0.0
        self.INERTIA = (1/12) * self.MASS * ((self.LENGTH)**2)
        self.DRAG_CONSTANT = 0.3
        self.CANARD_ONE_ORIENTATION = 0.0
        self.CANARD_TWO_ORIENTATION = 0.0
        self.CANARD_THREE_ORIENTATION = 0.0
        self.CANARD_FOUR_ORIENTATION = 0.0   
         
        
    def getInputs(self):
        """This takes in values for the wind speed in X & Y directions and the duration of the flight"""
        self.WSX = float(input("What is the windspeed in the X direction: "))
        self.WSY = float(input("What is the windspeed in the Y direction: "))
        self.DURATION = int(input("How many seconds would you like the rocket to fly for: "))
        return self
        
    def relativeWindSpeedCalc(self):
        self.RWS_X = self.WSX - self.VELOCITY_X
        self.RWS_Y = self.WSY - self.VELOCITY_Y
        return self.RWS_X, self.RWS_Y
        
    def torqueCalc(self): 
        """This method calculates the torque on the rocket
        
        may need to add a "self correcting" component to this"""
        #As the rockets orientation changes there will be less surface area in the direction of the wind to be hit
        crossSectionalArea_x = self.AREA * (np.cos(self.ORIENTATION_X)) 
        crossSectionalArea_y = self.AREA * (np.cos(self.ORIENTATION_Y))
        self.relativeWindSpeedCalc()
        force_x =  ((self.RWS_X)**2) * crossSectionalArea_x * self.DRAG_CONSTANT * np.sign(self.RWS_X)
        force_y =  ((self.RWS_Y)**2) * crossSectionalArea_y * self.DRAG_CONSTANT * np.sign(self.RWS_Y)
        torque_y =  self.DISTANCE * force_y
        torque_x =  self.DISTANCE * force_x
        return torque_x, torque_y

    def changeInAngularVelocity(self, torque_x, torque_y):
        """This method finds the change (in radians per second) of the rockets angular velocities"""
        self.ANGULAR_VELOCITY_X = self.TIMESTEP * torque_x / self.INERTIA
        self.ANGULAR_VELOCITY_Y = self.TIMESTEP * torque_y / self.INERTIA
        return self

    def changeInOrientation(self):
        """This method finds the change (in radians) of the rockets oreintation"""
        self.ORIENTATION_X += self.ANGULAR_VELOCITY_X * self.TIMESTEP
        self.ORIENTATION_Y += self.ANGULAR_VELOCITY_Y * self.TIMESTEP
        return self
    
    def changeInVelocity(self):
        #2 components - acceleration from rocket and acceleration from relative wind speed
        ACCEL_VERT = 30
        
        #THOUGHT:should these be swapped
        theta = self.ORIENTATION_X
        phi   = self.ORIENTATION_Y



        rotation_pitch = np.array([[1, 0, 0] ,
                                   [0, np.cos(theta), -1 * np.sin(theta)],
                                   [0, np.sin(theta), np.cos(theta)]])
        rotation_yaw   = np.array([[np.cos(phi), 0, np.sin(phi)],
                                   [0, 1, 0],
                                   [-1 * np.sin(phi), 0 , np.cos(phi)]])
        #This combines the rotation matrices so that
        #rotation_total = np.dot(rotation_pitch, rotation_yaw)
        rotation_total = np.dot(rotation_yaw,rotation_pitch)
        
        a_prime = np.array([0, 0, ACCEL_VERT])
        a = np.dot(rotation_total, a_prime)
 
        crossSectionalArea_X = self.AREA * (np.cos(self.ORIENTATION_X)) 
        crossSectionalArea_Y = self.AREA * (np.cos(self.ORIENTATION_Y))
        
        self.relativeWindSpeedCalc()
        force_X =  ((self.RWS_X)**2) * crossSectionalArea_X * self.DRAG_CONSTANT
        force_Y =  ((self.RWS_Y)**2) * crossSectionalArea_Y * self.DRAG_CONSTANT
        
        self.VELOCITY_X += (a[0] + force_X / self.MASS) * self.TIMESTEP
        self.VELOCITY_Y += (a[1] + force_Y / self.MASS) * self.TIMESTEP
        self.VELOCITY_Z += a[2] * self.TIMESTEP
        return self
    
    def changeInPosition(self):
        self.POSITION_X += self.VELOCITY_X * self.TIMESTEP
        self.POSITION_Y += self.VELOCITY_Y * self.TIMESTEP
        self.POSITION_Z += self.VELOCITY_Z * self.TIMESTEP
        return self
    
    
    def execute(self):
        time = 0
        Ps_X = []
        Ps_Y = []
        Ps_Z = []
        Os_X = []
        Os_Y = []
        Os_Z = []
        while (time < self.DURATION):
            Ps_X.append(self.POSITION_X)
            Ps_Y.append(self.POSITION_Y)
            Ps_Z.append(self.POSITION_Z)
            Os_X.append(self.ORIENTATION_X)
            Os_Y.append(self.ORIENTATION_Y)
            Os_Z.append(self.ORIENTATION_Z)
            
            self.changeInVelocity()
            self.changeInPosition()
            Tx, Ty = self.torqueCalc()
            self.changeInAngularVelocity(Tx, Ty)
            self.changeInOrientation()
            #if self.ORIENTATION_X > 0.3:
                #time = self.DURATION
            time += self.TIMESTEP
            
        return Ps_X, Ps_Y, Ps_Z, Os_X, Os_Y, Os_Z
        
""" OTHER FUNCTIONS - 1 for now"""

def plotter_3D(x,y,z):
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot3D(x,y,z, linewidth=5)
    ax.set_xlabel('X-axis') 
    ax.set_ylabel('Y-axis')  
    ax.set_zlabel('Z-axis')  
    plt.show()


"""MAIN"""        
    
rocket = Rocket()
rocket.getInputs()
Ps_X, Ps_Y, Ps_Z, Os_X, Os_Y, Os_Z = rocket.execute()
plotter_3D(Ps_X,Ps_Y,Ps_Z)
#plotter_3D(Os_X,Os_Y,Os_Z)

"""
THIS IS JUST MY THOUGHTS
Next steps:
consider that the force may be reduced as it is hitting at an angle (unsure)
consider the RWS in z direction (and different perpendicular distances in general?)



"""
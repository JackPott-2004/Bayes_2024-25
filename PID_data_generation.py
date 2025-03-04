#Rocket Sim with wind

"""
The rocket should never be more than ___ degrees from the vertical 

FINAL RESULT
you will be prompted how long you'd like the rocket to fly for
you will be prompted for the wind speed in x direction
you will be prompted for the wind speed in y direction
the rockets trajectory will be plotted, it will show it being cocked into the wind
optional - put a delay in when wind starts, variable wind speed with height

"""
 
#The imports
#from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation as R



class Rocket():
    
    def __init__(self):
        self.LENGTH = 1.0
        self.MASS = 5
        self.GRAVITY = -9.81
        self.AREA = 0.2
        self.TIMESTEP = 0.1
        #From centre of pressure to the centre of gravity,Should cock into positive wind values
        self.DISTANCE = 0.2 
        #initially vertical
        #Maybe I should put all of these into arrays
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
        self.INERTIA = (1/12) * self.MASS * ((self.LENGTH)**2) # abitrary but maybe the wrong formula
        self.DRAG_CONSTANT = 0.3
        self.CANARD_AREA = 1
        self.CANARD_DISTANCE = 0.2
        self.CANARDS_X_ORIENTATION = 0.0  # 
        self.CANARDS_Y_ORIENTATION = 0.0
        
    def getInputs(self):
        """This takes in values for the wind speed in X & Y directions and the duration of the flight"""
        #These wind speeds are what is coming at the rocket
        #self.WSX = -float(input("What is the windspeed in the X direction: "))
        #self.WSY = -float(input("What is the windspeed in the Y direction: "))
        #self.DURATION = int(input("How many seconds would you like the rocket to fly for: "))
        self.WSX = -10
        self.WSY = 0
        self.DURATION = 10
        return self
        
    def relativeWindSpeedCalc(self):
        self.RWS_X = self.VELOCITY_X - self.WSX  
        self.RWS_Y = self.VELOCITY_Y - self.WSY  
        self.RWS_Z = self.VELOCITY_Z
        return self
        
    def torqueCalc(self): 
        """This method calculates the torque on the rocket
        
        may need to add a "self correcting" component to this"""
        #As the rockets orientation changes there will be less surface area in the direction of the wind to be hit
        crossSectionalArea_x = self.AREA * (np.cos(self.ORIENTATION_X)) 
        crossSectionalArea_y = self.AREA * (np.cos(self.ORIENTATION_Y))
        self.relativeWindSpeedCalc()
        force_x =  ((self.RWS_X)**2) * crossSectionalArea_x * self.DRAG_CONSTANT * np.sign(self.RWS_X)
        force_y =  ((self.RWS_Y)**2) * crossSectionalArea_y * self.DRAG_CONSTANT * np.sign(self.RWS_Y)
        force_vector = np.array([force_x, force_y,0])
        moment_arm =np.array([0, 0, -self.DISTANCE])
        torque = np.cross(moment_arm, force_vector)
        print(torque)
        torque_x =  torque[0]
        torque_y =  torque[1]


        return torque_x, torque_y

    def canardTorqueCalc(self):
        #canards 1 and 3 & 2 and 4 are across from one another so will have the same area as they'll have the same orientation
        canards_X_area   = self.CANARD_AREA * np.sin(self.CANARDS_X_ORIENTATION + self.ORIENTATION_X) 
        canards_Y_area   = self.CANARD_AREA * np.sin(self.CANARDS_Y_ORIENTATION  + self.ORIENTATION_Y) 
        canardForce_x =  ((self.RWS_X)**2) * canards_X_area * self.DRAG_CONSTANT * np.sign(self.RWS_X)# this might wrong
        canardForce_y =  ((self.RWS_Y)**2) * canards_Y_area * self.DRAG_CONSTANT * np.sign(self.RWS_Y)
        canardForce_vector = np.array([canardForce_x, canardForce_y,0])
        moment_arm =np.array([0, 0, -self.DISTANCE])
        torque = np.cross(moment_arm, canardForce_vector)
        canardTorque_x =  torque[0]
        canardTorque_y =  torque[1]
        return canardTorque_x, canardTorque_y

    def canardTorqueCalc_WithSetCanards(self,X_orientation, Y_orientation ):
        canards_X_area   = self.CANARD_AREA * np.sin(self.CANARDS_X_ORIENTATION  + self.ORIENTATION_X) * np.sin(self.ORIENTATION_Y)
        canards_Y_area   = self.CANARD_AREA * np.sin(self.CANARDS_Y_ORIENTATION  + self.ORIENTATION_Y) * np.sin(self.ORIENTATION_X)
        canardForce_x =  ((self.RWS_X)**2) * canards_X_area * self.DRAG_CONSTANT * np.sign(self.RWS_X)# this might wrong
        canardForce_y =  ((self.RWS_Y)**2) * canards_Y_area * self.DRAG_CONSTANT * np.sign(self.RWS_Y)
        canardForce_vector = np.array([canardForce_x, canardForce_y,0])
        moment_arm =np.array([0, 0, -self.DISTANCE])
        torque = np.cross(moment_arm, canardForce_vector)
        canardTorque_x =  torque[0]
        canardTorque_y =  torque[1]
        return canardTorque_x, canardTorque_y

    def changeInAngularVelocity(self, torque_x, torque_y, canardTorque_x, canardTorque_y):
        """This method finds the change (in radians per second) of the rockets angular velocities"""
        self.ANGULAR_VELOCITY_X = self.TIMESTEP * (torque_x + canardTorque_x) / self.INERTIA # here too
        self.ANGULAR_VELOCITY_Y = self.TIMESTEP * (torque_y + canardTorque_y)/ self.INERTIA
        return self

    def changeInOrientation(self):
        """This method finds the change (in radians) of the rockets oreintation"""
        self.ORIENTATION_X += self.ANGULAR_VELOCITY_X * self.TIMESTEP
        self.ORIENTATION_Y += self.ANGULAR_VELOCITY_Y * self.TIMESTEP
        return self
    
    #Is this method of numerical integration alright?
    #What method is it
    #euler,euler - cromer?
    def changeInVelocity(self):
        #2 components - acceleration from rocket and acceleration from relative wind speed
        ACCEL_VERT = 30
        accel_array = [0,0,ACCEL_VERT]

        theta = self.ORIENTATION_X
        phi   = self.ORIENTATION_Y

        axis_x = np.array([1,0,0])
        axis_y = np.array([0,1,0])

        q_x = R.from_rotvec(theta * axis_x)
        q_y = R.from_rotvec(phi * axis_y)

        total_rotation = q_x * q_y
        a = total_rotation.apply(accel_array)

        crossSectionalArea_X = self.AREA * (np.cos(self.ORIENTATION_X)) 
        crossSectionalArea_Y = self.AREA * (np.cos(self.ORIENTATION_Y))
        
        self.relativeWindSpeedCalc()
        force_X =  ((self.RWS_X)**2) * crossSectionalArea_X * self.DRAG_CONSTANT
        force_Y =  ((self.RWS_Y)**2) * crossSectionalArea_Y * self.DRAG_CONSTANT
        
        force_array = np.array([force_X, force_Y,0])
        force_rotated = total_rotation.apply(force_array)

        self.VELOCITY_X += (a[0] + force_rotated[0] / self.MASS) * self.TIMESTEP
        self.VELOCITY_Y += (a[1] + force_rotated[1] / self.MASS) * self.TIMESTEP
        self.VELOCITY_Z += (a[2] + self.GRAVITY) * self.TIMESTEP
        return self
    
    def changeInPosition(self):
        self.POSITION_X += self.VELOCITY_X * self.TIMESTEP
        self.POSITION_Y += self.VELOCITY_Y * self.TIMESTEP
        self.POSITION_Z += self.VELOCITY_Z * self.TIMESTEP
        return self
    
    def Proportional(self,orientation_x, orientation_y):
        if orientation_x == 0 and orientation_y == 0:
            return 0, 0
        else:
            NUM_ANGLES = 120
            X_trackerArray = np.empty(NUM_ANGLES)
            Y_trackerArray = np.empty(NUM_ANGLES)
            ANGLES_ARRAY = np.empty(NUM_ANGLES) 
            TORQUES_X = np.empty(NUM_ANGLES)
            TORQUES_Y = np.empty(NUM_ANGLES)

            for i in range(NUM_ANGLES):
                ANGLES_ARRAY[i] = (-np.pi + np.pi * 2 * i / NUM_ANGLES)
                X_trackerArray[i] = i
                Y_trackerArray[i] = i
            for j in range(NUM_ANGLES):
                tx, ty = self.canardTorqueCalc_WithSetCanards(ANGLES_ARRAY[j],ANGLES_ARRAY[j])
                TORQUES_X[j] = tx
                TORQUES_Y[j] = ty

            if orientation_x < 0 :# may use np.sign here to reduce no. lines
                for m in range(NUM_ANGLES):
                    if TORQUES_X[m] > 0:
                        np.delete(TORQUES_X[m],m)
                        np.delete(X_trackerArray[m],m)
            if orientation_x > 0 :# may use np.sign here to reduce no. lines
                for n in range(NUM_ANGLES):
                    if TORQUES_X[n] < 0:
                        np.delete(TORQUES_X[n],n)
                        np.delete(X_trackerArray[n],n)
            if orientation_y < 0 :# may use np.sign here to reduce no. lines
                for o in range(NUM_ANGLES):
                    if TORQUES_X[o] < 0:
                        np.delete(TORQUES_Y[o],o)
                        np.delete(Y_trackerArray[o],o)
            if orientation_y > 0 :# may use np.sign here to reduce no. lines
                for p in range(NUM_ANGLES):
                    if TORQUES_X[p] > 0:
                        np.delete(TORQUES_Y[p],p)
                        np.delete(Y_trackerArray[p],p)

            X_TORQUES_SORTED = np.sort(TORQUES_X)
            Y_TORQUES_SORTED = np.sort(TORQUES_Y)     
            round_to = round(0.8 * NUM_ANGLES) #Kameran said to take the "20th" best one is there was 100
            return X_TORQUES_SORTED[-1], Y_TORQUES_SORTED[-1]

    def error(self): #If i wanted to add a variable gain term - might not as if it works well it would be unneccessary
        DESIRED_AXIS = [0,0,1]
        AXIS = []

        return
    
    def execute(self):
        time = 0
        Ps_X = []
        Ps_Y = []
        Ps_Z = []
        Os_X = []
        Os_Y = []
        Os_Z = [] 
        xCanardAgles = []
        yCanardAgles = []

        while (time < self.DURATION):
            Ps_X.append(self.POSITION_X)
            Ps_Y.append(self.POSITION_Y)
            Ps_Z.append(self.POSITION_Z)
            Os_X.append(self.ORIENTATION_X)
            Os_Y.append(self.ORIENTATION_Y) 
            Os_Z.append(self.ORIENTATION_Z) 
            xCanardAgles.append(self.CANARDS_X_ORIENTATION)
            yCanardAgles.append(self.CANARDS_Y_ORIENTATION)
            self.changeInVelocity()
            self.changeInPosition()
            Tx, Ty = self.torqueCalc()
            CTx, CTy = self.canardTorqueCalc()
            self.changeInAngularVelocity(Tx, Ty, CTx, CTy)
            self.changeInOrientation()
            self.CANARDS_X_ORIENTATION, self.CANARDS_Y_ORIENTATION = self.Proportional(self.ORIENTATION_X, self.ORIENTATION_Y)
            #if self.ORIENTATION_X > 0.3: there also would be one one of these IFs for the y axis, as we aren't allowed to control the rocket when it is out of control
                #time = self.DURATION
            time += self.TIMESTEP
            
        return Ps_X, Ps_Y, Ps_Z, Os_X, Os_Y, Os_Z, xCanardAgles,yCanardAgles
        
""" FUNCTIONS """

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

"""MAIN"""        
def main():    
    rocket = Rocket()
    rocket.getInputs()
    times = np.linspace(0,rocket.DURATION,int(rocket.DURATION/rocket.TIMESTEP)+1)
    Ps_X, Ps_Y, Ps_Z, Os_X, Os_Y, Os_Z, XCA, YCA = rocket.execute()
    plotter_3D(Ps_X,Ps_Y,Ps_Z)
    #plotter_3D(Os_X,Os_Y,Os_Z)
    #plot_2d(XCA,times) #Not sure why it isnt showing the plot
    #plot_2d(Os_Y,times)
main()

"""
BIG STEPS TO TAKE

experiment with inheritance


"""

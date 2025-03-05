#Rocket Sim with wind

"""
The rocket should never be more than ___ degrees from the vertical 

FINAL RESULT
you will be prompted how long you'd like the rocket to fly for
you will be prompted for the wind speed in x direction
you will be prompted for the wind speed in y direction
the rockets trajectory will be plotted, it will show it being cocked into the wind
but when the proportional canards are implemented it will just drift because of the wind and no cocking
optional - put a delay in when wind starts, variable wind speed with height

"""
#from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation as R



class Rocket():
    
    def __init__(self):
        #The intrinsics
        self.LENGTH = 1.0
        self.MASS = 5
        self.GRAVITY = -9.81
        self.AREA = 0.2
        self.TIMESTEP = 0.1
        self.DISTANCE = 0.2
        self.INERTIA = 0.1
        self.DRAG_CONSTANT = 0.3
        self.CANARD_AREA = 1
        self.CANARD_DISTANCE = 0.2 
        #variable ones
        self.ORIENTATIONS = np.zeros(3)
        self.ANGULAR_VELOCITies = np.zeros(3)
        self.VELOCITIES = np.zeros(3)
        self.POSITIONS = np.zeros(3)
       
        self.CANARDS_X_ORIENTATION = 0.0  # 
        self.CANARDS_Y_ORIENTATION = 0.0
        
    def getInputs(self):
        """This takes in values for the wind speed in X & Y directions and the duration of the flight"""
        #These wind speeds are what is coming at the rocket
        #self.WSX = -float(input("What is the windspeed in the X direction: "))
        #self.WSY = -float(input("What is the windspeed in the Y direction: "))
        #self.DURATION = int(input("How many seconds would you like the rocket to fly for: "))
        self.WSX = -10 #For ease of testing
        self.WSY = 0
        self.DURATION = 10
        return self
        
    def relativeWindSpeedCalc(self):
        windSpeed = [self.WSX,self.WSY,0]
        self.RWS = - self.VELOCITIES - windSpeed
        return self
        
    def torqueCalc(self): 
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
    """
    def canardTorqueCalc(self):
        #canards 1 and 3 & 2 and 4 are across from one another so will have the same area as they'll have the same orientation
        canards_X_area   = self.CANARD_AREA * np.sin(-self.CANARDS_X_ORIENTATION + self.ORIENTATIONS[0]) 
        canards_Y_area   = self.CANARD_AREA * np.sin(-self.CANARDS_Y_ORIENTATION  + self.ORIENTATIONS[1]) 
        canardForce_x =  ((self.RWS[0])**2) * canards_X_area * self.DRAG_CONSTANT * np.sign(self.RWS[0])# this might wrong
        canardForce_y =  ((self.RWS[1])**2) * canards_Y_area * self.DRAG_CONSTANT * np.sign(self.RWS[1])
        canardForce_vector = np.array([canardForce_x, canardForce_y,0])
        moment_arm =np.array([0, 0, -self.DISTANCE])
        torque = np.cross(moment_arm, canardForce_vector)
        canardTorque_x =  torque[0]
        canardTorque_y =  torque[1]
        return canardTorque_x, canardTorque_y"""
        
    """
    def canardTorqueCalc_WithSetCanards(self,X_orientation, Y_orientation ):#i dislike that these are seperate^^^^
        canards_X_area   = self.CANARD_AREA * np.sin(X_orientation  + self.ORIENTATIONS[0]) * np.sin(self.ORIENTATIONS[1])
        canards_Y_area   = self.CANARD_AREA * np.sin(Y_orientation  + self.ORIENTATIONS[1]) * np.sin(self.ORIENTATIONS[0])
        canardForce_x =  ((self.RWS[0])**2) * canards_X_area * self.DRAG_CONSTANT * np.sign(self.RWS[0])# this might wrong
        canardForce_y =  ((self.RWS[1])**2) * canards_Y_area * self.DRAG_CONSTANT * np.sign(self.RWS[1])
        canardForce_vector = np.array([canardForce_x, canardForce_y,0])
        moment_arm =np.array([0, 0, -self.DISTANCE])
        torque = np.cross(moment_arm, canardForce_vector)
        canardTorque_x =  torque[0]
        canardTorque_y =  torque[1]
        return canardTorque_x, canardTorque_y
    """

    def changeInAngularVelocity(self, torque_x, torque_y):
        """This method finds the change (in radians per second) of the rockets angular velocities"""
        self.ANGULAR_VELOCITY_X = self.TIMESTEP * (torque_x) / self.INERTIA # here too
        self.ANGULAR_VELOCITY_Y = self.TIMESTEP * (torque_y)/ self.INERTIA
        return self

    def changeInOrientation(self):
        """This method finds the change (in radians) of the rockets oreintation"""
        self.ORIENTATIONS[0] += self.ANGULAR_VELOCITY_X * self.TIMESTEP
        self.ORIENTATIONS[1] += self.ANGULAR_VELOCITY_Y * self.TIMESTEP
        return self
    
    #Is this method of numerical integration alright?
    #What method is it
    #euler,euler - cromer?
    def changeInVelocity(self):
        ACCEL_VERT = 30
        accel_array = [0,0,ACCEL_VERT]

        theta = self.ORIENTATIONS[0]
        phi   = self.ORIENTATIONS[1]

        axis_x = np.array([1,0,0])
        axis_y = np.array([0,1,0])

        q_x = R.from_rotvec(theta * axis_x)
        q_y = R.from_rotvec(phi * axis_y)

        total_rotation = q_x * q_y
        a = total_rotation.apply(accel_array)

        crossSectionalArea_X = self.AREA * (np.cos(theta)) 
        crossSectionalArea_Y = self.AREA * (np.cos(phi))
        
        
        self.relativeWindSpeedCalc()
        force_X =  ((self.RWS[0])**2) * crossSectionalArea_X * self.DRAG_CONSTANT
        force_Y =  ((self.RWS[1])**2) * crossSectionalArea_Y * self.DRAG_CONSTANT
        
        force_array = np.array([force_X, force_Y,0])
        force_rotated = total_rotation.apply(force_array)

        self.VELOCITIES[0] += (a[0] + force_rotated[0] / self.MASS) * self.TIMESTEP
        self.VELOCITIES[1] += (a[1] + force_rotated[1] / self.MASS) * self.TIMESTEP
        self.VELOCITIES[2] += (a[2] + self.GRAVITY) * self.TIMESTEP
        return self
    
    def changeInPosition(self):
        self.POSITIONS[0] += self.VELOCITIES[0] * self.TIMESTEP
        self.POSITIONS[1] += self.VELOCITIES[1] * self.TIMESTEP
        self.POSITIONS[2] += self.VELOCITIES[2] * self.TIMESTEP
        return self
    """
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
            #Dont think below considers enough rotations
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
    """

    """
    def error(self): #If i wanted to add a variable gain term - might not as if it works well it would be unneccessary
        DESIRED_AXIS = [0,0,1]
        AXIS = []

        return
    """
    
    def execute(self):
        time = 0
        Ps = [[],[],[]]
        Os = [[],[],[]] 
        xCanardAgles = []
        yCanardAgles = []

        while (time < self.DURATION):

            Ps[0].append(self.POSITIONS[0])
            Ps[1].append(self.POSITIONS[1])
            Ps[2].append(self.POSITIONS[2])
            Os[0].append(self.ORIENTATIONS[0])
            Os[1].append(self.ORIENTATIONS[1]) 
            Os[2].append(self.ORIENTATIONS[2]) 
            xCanardAgles.append(self.CANARDS_X_ORIENTATION)
            yCanardAgles.append(self.CANARDS_Y_ORIENTATION)
            self.changeInVelocity()
            self.changeInPosition()
            Tx, Ty = self.torqueCalc()
            #CTx, CTy = self.canardTorqueCalc()
            self.changeInAngularVelocity(Tx, Ty)
            self.changeInOrientation()
            #self.CANARDS_X_ORIENTATION, self.CANARDS_Y_ORIENTATION = self.Proportional(self.ORIENTATIONS[0], self.ORIENTATIONS[1])
            #if self.ORIENTATION_X > 0.3: there also would be one one of these IFs for the y axis, as we aren't allowed to control the rocket when it is out of control
                #time = self.DURATION
            time += self.TIMESTEP
            
        return Ps, Os, xCanardAgles, yCanardAgles
        
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
    Ps, Os, XCA, YCA = rocket.execute()
    plotter_3D(Ps[0],Ps[1],Ps[2])
    #plotter_3D(Os[0],Os[1],Os[2])
    #plot_2d(XCA,times) #Not sure why it isnt showing the plot
    #plot_2d(Os[1],times)
main()

"""
BIG STEPS TO TAKE
migrate variables to arrays

experiment with inheritance


"""

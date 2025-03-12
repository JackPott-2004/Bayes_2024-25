
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation as R
import math

def rotateToRocketsAxis(array):
        theta = np.pi/2
        phi = np.pi/2
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

accel_array = [0, 0, 30]
rotated = rotateToRocketsAxis(accel_array)
print(rotated)
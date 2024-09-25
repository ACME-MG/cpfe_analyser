"""
 Title:         Texture
 Description:   Analyses the texture of different meshes
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
from scipy.special import sph_harm

def euler_to_spherical(phi1, Phi, phi2):
    """
    Convert Euler-Bunge angles to spherical coordinates.
    
    Parameters:
        phi1 (float): First Euler angle (rotation around Z-axis).
        Phi (float): Second Euler angle (rotation around X-axis).
        phi2 (float): Third Euler angle (rotation around new Z'-axis).
        
    Returns:
        (theta, phi): Colatitudinal and azimuthal angles in radians.
    """
    # Rotation matrices
    R_z1 = np.array([[np.cos(phi1), -np.sin(phi1), 0], [np.sin(phi1),  np.cos(phi1), 0], [0, 0, 1]])
    R_x = np.array([[1, 0, 0], [0, np.cos(Phi), -np.sin(Phi)], [0, np.sin(Phi), np.cos(Phi)]])
    R_z2 = np.array([[np.cos(phi2), -np.sin(phi2), 0], [np.sin(phi2),  np.cos(phi2), 0],  [0, 0, 1]])
    R = R_z2 @ R_x @ R_z1
    
    # Transform the unit vector [1, 0, 0]
    transformed_vector = R @ np.array([1, 0, 0])
    
    # Compute spherical coordinates
    x, y, z = transformed_vector
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)  # Colatitudinal angle
    phi = np.arctan2(y, x)    # Azimuthal angle
    
    return theta, phi

def get_texture_index(orientations:list, weights:list) -> list:
    """
    Calculates the orientation distribution function
    
    Parameters:
    * `orientations`: List of orientations in euler-bunge (rads)
    * `weights`:      List of weights (e.g., grain volume)
    
    Returns the ODF
    """
    n_max = 4
    odf = np.zeros((n_max+1, 2*n_max+1), dtype=complex)
    for orientation, weight in zip(orientations, weights):
        theta, phi = euler_to_spherical(*orientation)
        theta = orientation[0] # azimuthal (longitudinal) angle
        phi   = orientation[1] # polar (colatitudinal) angle
        for n in range(n_max+1):
            for m in range(-n, n+1):
                Y_lm = sph_harm(m, n, theta, phi)
                odf[n, m+n] += weight*Y_lm
    total_weight = np.sum(weights)
    odf /= total_weight
    texture_index = 0.0
    for n in range(n_max+1):
        for m in range(-n, n+1):
            texture_index += np.abs(odf[n, m+n])**2
    return texture_index

# orientations = [[0, np.pi/3, np.pi/4], [np.pi/4, 2*np.pi/3, np.pi/4]]
# weights = [1, 1]
# # texture_index = get_texture_index(orientations, weights)
orientations = [(np.pi/3, np.pi/4), (2*np.pi/3, np.pi/4)]  # Example spherical coordinates
weights = [1, 5]
texture_index = get_texture_index(orientations, weights)
print("Texture Index:", texture_index)
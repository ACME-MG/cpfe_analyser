"""
 Title:         Manual
 Description:   Manual calculation of elastic strain
 References:    [1] https://www.researchgate.net/publication/324088567_Computing_Euler_angles_with_Bunge_convention_from_rotation_matrix
 Author:        Janzen Choi

"""

# Libraries
import itertools, math, numpy as np
import matplotlib.pyplot as plt
from neml.math import rotations, tensors
from neml import elasticity

def reorient_euler(euler:list) -> list:
    """
    Inverts the euler angle from passive/active to active/passive

    Parameters:
    * `euler`: The euler angle (rads)

    Returns the inverted euler angle
    """
    orientation = rotations.CrystalOrientation(euler[0], euler[1], euler[2], angle_type="radians", convention="bunge")
    inverse = orientation.inverse()
    new_euler = inverse.to_euler(angle_type="radians", convention="bunge")
    return list(new_euler)

def multiply_tensors(tensor_1:np.array, tensor_2:np.array) -> np.array:
    """
    Multiplies a 3x3x3x3 tensor by a 3x3 tensor
    
    Parameters:
    * `tensor_1`: The 3x3 tensor
    * `tensor_2`: The 3x3x3x3 tensor
    
    Returns the product tensor as a 3x3 array
    """
    tensor_3 = np.zeros((3, 3))    
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    tensor_3[i,j] += tensor_2[i,j,k,l]*tensor_1[k,l]
    return tensor_3

def deg_to_rad(degrees:float) -> float:
    """
    Converts degrees to radians

    Parameters:
    * `degrees`: The degrees to be converted

    Returns the converted radians
    """
    if isinstance(degrees, list):
        return [deg_to_rad(d) for d in degrees]
    return degrees * math.pi / 180

def sym_tensor_part(C):
    """
    Take a Mandel stiffness in my notation and convert it back to a full tensor
    """
    mandel = ((0,0),(1,1),(2,2),(1,2),(0,2),(0,1))
    mandel_mults = (1,1,1,np.sqrt(2),np.sqrt(2),np.sqrt(2))
    Ct = np.zeros((3,3,3,3))
    for a in range(6):
        for b in range(6):
            ind_a = itertools.permutations(mandel[a], r=2)
            ind_b = itertools.permutations(mandel[b], r=2)
            ma = mandel_mults[a]
            mb = mandel_mults[b]
            indexes = tuple(ai+bi for ai, bi in itertools.product(ind_a, ind_b))
            for ind in indexes:
                Ct[ind] = C[a,b] / ma*mb
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    if l < k:
                        Ct[i,j,k,l] = 0.0
    return Ct

def orientation_to_c(orientation:list) -> np.array:
    """
    Converts an orientation to the elasticity tensor
    in the reference frame
    
    Parameters:
    * `orientation`: The passive orientation
    
    Returns the elasticity tensor
    """
    euler = deg_to_rad(orientation)
    active_euler = reorient_euler(euler) # active
    active_orientation = rotations.CrystalOrientation(*active_euler, angle_type="radians", convention="bunge")
    e_model = elasticity.CubicLinearElasticModel(205000, 138000, 126000, "components")
    c_6x6 = e_model.C_tensor(20, active_orientation)
    c_6x6 = c_6x6.data.reshape(6,6)
    c_3x3x3x3 = sym_tensor_part(c_6x6)
    return c_3x3x3x3

# Main function
if __name__ == "__main__":

    # Define orientations
    ORIENTATIONS = [[0,135,45], [68,38,72], [110,75,80], [0,0,0]]
    COLOUR_LIST = ["green", "black", "blue", "red"]
    LABEL_LIST = ["220", "111", "311", "200"]

    # Format plot
    plt.figure(figsize=(5,5))
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    plt.xlabel("Elastic Strain (μmm/mm)")
    plt.ylabel("Applied Stress (MPa)")

    # Iterate through grains
    stress_list = list(np.linspace(0,1000,100))
    for i, orientation in enumerate(ORIENTATIONS):
        c = orientation_to_c(orientation)
        strain_list = []
        for stress in stress_list:
            stress_tensor = np.array([[stress,0,0], [0,0,0], [0,0,0]])
            strain_tensor = multiply_tensors(stress_tensor, c)
            strain_list.append(strain_tensor[0,0]*1e6)
        plt.plot(strain_list, stress_list, color=COLOUR_LIST[i], label=LABEL_LIST[i])

    # Save the plot
    plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white")
    plt.savefig("plot_es.png")


# Libraries
import itertools, math, numpy as np
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

def get_c(c_11:float, c_12:float, c_44:float) -> np.array:
    """
    Gets the elasticity tensor
    
    Parameters:
    * `c_11`: The 11 element in the elasticity tensor
    * `c_12`: The 12 element in the elasticity tensor
    * `c_44`: The 44 element in the elasticity tensor
    
    Returns the elasticity tensor as a 3x3x3x3 array
    """
    c = np.zeros((3, 3, 3, 3))
    for i in range(3):
        c[i, i, i, i] = c_11
        for j in range(3):
            if i != j:
                c[i, i, j, j] = c_12
                c[i, j, i, j] = c_44
    return c

def euler_to_matrix(euler:list) -> np.array:
    """
    Determines the orientation matrix of a set of euler-bunge angles (rads);
    from Ref. [1]
    
    Parameters:
    * `euler`: The euler angle in euler-bunge form

    Returns the 3x3 orientation matrix as an array
    """
    om_11 = math.cos(euler[0])*math.cos(euler[2]) - math.sin(euler[0])*math.sin(euler[2])*math.cos(euler[1])
    om_12 = math.sin(euler[0])*math.cos(euler[2]) + math.cos(euler[0])*math.sin(euler[2])*math.cos(euler[1])
    om_13 = math.sin(euler[2])*math.sin(euler[1])
    om_21 = -math.cos(euler[0])*math.sin(euler[2]) - math.sin(euler[0])*math.cos(euler[2])*math.cos(euler[1])
    om_22 = -math.sin(euler[0])*math.sin(euler[2]) + math.cos(euler[0])*math.cos(euler[2])*math.cos(euler[1])
    om_23 = math.cos(euler[2])*math.sin(euler[1])
    om_31 = math.sin(euler[0])*math.sin(euler[1])
    om_32 = -math.cos(euler[0])*math.sin(euler[1])
    om_33 = math.cos(euler[1])
    om = np.array([[om_11, om_12, om_13], [om_21, om_22, om_23], [om_31, om_32, om_33]])
    return om

def get_c_prime(r:np.array, c:np.array) -> list:
    """
    Rotates the elasticity tensor
    
    Parameters:
    * `r`: The 3x3 rotation matrix
    * `c`: The 3x3x3x3 elasticity tensor
    
    Returns the rotated elasticity tensor as a 3x3x3x3 array
    """
    c_prime = np.zeros((3, 3, 3, 3))    
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            for o in range(3):
                                for p in range(3):
                                    r_factor = r[i][m]*r[j][n]*r[k][o]*r[l][p]
                                    c_prime[i,j,k,l] += r_factor*c[m,n,o,p]
    return c_prime

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

def invert_3x3x3x3(tensor:np.array) -> np.array:
    """
    Inverts a 3x3x3x3 tensor
    
    Parameters:
    * `tensor`: The tensor to be inverted
    
    Returns the inverted tensor
    """
    tensor_inv = np.empty_like(tensor)
    for i in range(tensor.shape[0]):
        for j in range(tensor.shape[1]):
            try:
                tensor_inv[i,j] = np.linalg.inv(tensor[i,j])
            except np.linalg.LinAlgError:
                tensor_inv[i,j] = np.nan
    return tensor_inv

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

def orientation_to_c(euler:list) -> np.array:
    """
    Converts an orientation to the elasticity tensor
    in the reference frame, inverted
    
    Parameters:
    * `euler`: The passive orientation (degrees)
    
    Returns the elasticity tensor
    """
    euler = deg_to_rad(euler)
    euler = reorient_euler(euler) # active
    r = euler_to_matrix(euler)
    c = get_c(205000, 138000, 126000) # MPa
    c_p = get_c_prime(r, c)
    c_pi = invert_3x3x3x3(c_p)
    return c_pi

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

# Define the orientation
euler = [0,135,45]

# Calculate C tensor manually
active_euler = reorient_euler(deg_to_rad(euler)) # active
r = euler_to_matrix(active_euler)
c = get_c(205000, 138000, 126000) # MPa
c_p = get_c_prime(r, c)
c_inverse_1 = invert_3x3x3x3(c_p)

# Calculate C tensor using NEML code
e_model = elasticity.CubicLinearElasticModel(205000, 138000, 126000, "components")
passive_orientation = rotations.Orientation(*euler, angle_type="degrees", convention="bunge")
active_orientation = passive_orientation.inverse()
c_inverse = e_model.C_tensor(20, active_orientation).inverse()
c_inverse = c_inverse.data.reshape(6, 6)
c_inverse_2 = sym_tensor_part(c_inverse)

print(c_inverse_1[0,0])
print("=======================================================")
print(c_inverse_2[0,0])

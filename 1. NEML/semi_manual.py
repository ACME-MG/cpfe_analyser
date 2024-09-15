
# Libraries
import itertools, math, numpy as np
import matplotlib.pyplot as plt
from neml.cp import polycrystal, crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polefigures
from neml.math import rotations, nemlmath, tensors
from neml import elasticity

def transpose(list_of_lists:list) -> list:
    """
    Transposes a 2D list of lists
    
    Parameters:
    * `list_of_lists`: A list of lists (i.e., a 2D grid)
    
    Returns the transposed list of lists
    """
    transposed = np.array(list_of_lists).T.tolist()
    return transposed

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

# Model parameters
TAU_SAT  = 108.35
B        = 0.5840
TAU_0    = 120.21*1e6 # CRSS is super high
GAMMA_0  = 3.3333e-5
N        = 2.5832

# Define orientations
ORIENTATIONS = [[0,135,45], [68,38,72], [110,75,80], [0,0,0]]
COLOUR_LIST  = ["green", "black", "blue", "red"]
LABEL_LIST   = ["220", "111", "311", "200"]

# Define orientations
passive_eulers = deg_to_rad(ORIENTATIONS) # passive
active_eulers = [reorient_euler(euler) for euler in passive_eulers] # active
orientations = [rotations.CrystalOrientation(*euler, angle_type="radians", convention="bunge") for euler in active_eulers]
weights = [1] * len(orientations)

# Define lattice
lattice = crystallography.CubicLattice(1.0)
lattice.add_slip_system([1,1,0], [1,1,1])

# Define model
e_model    = elasticity.CubicLinearElasticModel(284039, 121731, 162000, "components")
# e_model    = elasticity.CubicLinearElasticModel(205000, 138000, 126000, "components")
str_model  = slipharden.VoceSlipHardening(TAU_SAT, B, TAU_0)
slip_model = sliprules.PowerLawSlipRule(str_model, GAMMA_0, N)
i_model    = inelasticity.AsaroInelasticity(slip_model)
k_model    = kinematics.StandardKinematicModel(e_model, i_model)
sc_model   = singlecrystal.SingleCrystalModel(k_model, lattice, miter=16, max_divide=2, verbose=False)
pc_model   = polycrystal.TaylorModel(sc_model, orientations, nthreads=4, weights=weights)

# Define driver control parameters
strain_tensor = np.array([[1.0,0,0], [0,-0.5,0], [0,0,-0.5]])
strain_rate = 1.0e-4
strain_tensor *= strain_rate
steps = 500
strain_max = 0.01
temperature = 20
dt = strain_max / steps / strain_rate

# Define driver auxiliary parameters
h_n = pc_model.init_store()
d_inc = nemlmath.sym(0.5*(strain_tensor+strain_tensor.T))
w_inc = nemlmath.skew(0.5*(strain_tensor-strain_tensor.T))
s_n = np.zeros((6,))
d_n = np.zeros((6,))
w_n = np.zeros((3,))
u_n = 0.0
p_n = 0.0

# Initialise storage
results = {"history": [], "stress": []}

# Run custom driver
for i in range(steps):
    
    # Run step
    d_np1 = d_n + d_inc * dt
    w_np1 = w_n + w_inc * dt
    s_np1, h_np1, A_np1, B_np1, u_np1, p_np1 = pc_model.update_ld_inc(d_np1, d_n, w_np1, w_n, temperature,
                                                                      temperature, dt, 0, s_n, h_n, u_n, p_n)
    
    # Update variables
    d_n = np.copy(d_np1)
    w_n = np.copy(w_np1)
    s_n = np.copy(s_np1)
    h_n = np.copy(h_np1)
    u_n = u_np1
    p_n = p_np1
    
    # Store results
    results["history"].append(h_n)
    results["stress"].append(s_n)

# Extract results
history = np.array(results["history"])
stress_list = np.array(results["stress"])[:,0]
nhist = sc_model.nstore

# Read stress tensors
stress_xx = history[:, pc_model.n*nhist+0:pc_model.n*nhist+6*pc_model.n:6]
stress_yy = history[:, pc_model.n*nhist+1:pc_model.n*nhist+6*pc_model.n:6]
stress_zz = history[:, pc_model.n*nhist+2:pc_model.n*nhist+6*pc_model.n:6]
stress_xy = history[:, pc_model.n*nhist+3:pc_model.n*nhist+6*pc_model.n:6]
stress_yz = history[:, pc_model.n*nhist+4:pc_model.n*nhist+6*pc_model.n:6]
stress_xz = history[:, pc_model.n*nhist+5:pc_model.n*nhist+6*pc_model.n:6]

# Create compiled stress tensor
stress_tensor = np.zeros((len(stress_xx), len(stress_xx[0]), 3, 3))
for i in range(len(stress_xx)):
    for j in range(len(stress_xx[0])):
        a = np.sqrt(2)
        stress_tensor[i,j] = np.array([
            [stress_xx[i,j], stress_xy[i,j]/a, stress_xz[i,j]]/a,
            [stress_xy[i,j]/a, stress_yy[i,j], stress_yz[i,j]]/a,
            [stress_xz[i,j]/a, stress_yz[i,j]/a, stress_zz[i,j]]
        ])

# Get elastic strain tensors
elastic_tensor = np.zeros((len(stress_xx), len(stress_xx[0]), 6))
for i, state in enumerate(history):
    orientation_list = pc_model.orientations(state)
    for j, orientation in enumerate(orientation_list):
        
        # new_e_model = elasticity.CubicLinearElasticModel(284039, 121731, 162000, "components")
        new_e_model = elasticity.CubicLinearElasticModel(205000, 138000, 126000, "components")
        
        # c_inverse_6x6 = new_e_model.C_tensor(20, orientation.inverse()).inverse()
        # c_inverse_6x6 = c_inverse_6x6.data.reshape(6,6)
        # c_inverse_3x3x3x3 = sym_tensor_part(c_inverse_6x6)
        # strain_tensor = multiply_tensors(stress_tensor[i,j], c_inverse_3x3x3x3)
        # strain_tensor = tensors.Symmetric(strain_tensor).data
        # elastic_tensor[i,j] = strain_tensor
        
        tensor = new_e_model.S_tensor(20, orientation.inverse()) # active
        strain_tensor = tensor * tensors.Symmetric(stress_tensor[i,j])
        elastic_tensor[i,j] = strain_tensor.data

elastic_xx = elastic_tensor[:,:,0]
elastic_xx_history = transpose(elastic_xx)

# Initialise plot
plt.figure(figsize=(5,5))
plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":")
plt.xlabel("Elastic Strain (μmm/mm)")
plt.ylabel("Applied Stress (MPa)")

# Iterate through grains
for i in range(len(ORIENTATIONS)):
    elastic_strain = transpose(elastic_xx_history[i])
    elastic_strain = [1e6*es for es in elastic_strain]
    plt.plot(elastic_strain, stress_list, color=COLOUR_LIST[i], label=LABEL_LIST[i])

# Save the plot
plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white")
plt.savefig("plot_es.png")

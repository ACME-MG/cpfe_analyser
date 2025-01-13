"""
 Title:         Manual
 Description:   Manual calculation of elastic strain
 References:    [1] https://www.researchgate.net/publication/324088567_Computing_Euler_angles_with_Bunge_convention_from_rotation_matrix
 Author:        Janzen Choi

"""

# Libraries
import itertools
import math, numpy as np
import matplotlib.pyplot as plt
from neml import elasticity, drivers
from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal
from neml.math import rotations, tensors

# Define families
COLOUR_LIST = ["green", "black", "blue", "red"]
FAMILY_LIST = [[2,2,0], [1,1,1], [3,1,1], [2,0,0]]

# Constants
NUM_THREADS = 8
STRAIN_RATE = 1e-4
MAX_STRAIN  = 0.1

# Model parameters
TAU_SAT  = 825
B        = 2 # 0.3
TAU_0    = 112
GAMMA_0  = 3.3333e-5
N        = 15

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

def get_grain_family(orientations:list, crystal_direction:list,
                     sample_direction:list, threshold:float=10.0) -> list:
    """
    Groups a list of orientations to a family

    Parameters:
    * `orientations`:      The list of passive euler-bunge angles (rads)
    * `crystal_direction`: The crystal direction
    * `sample_direction`:  The sample direction
    * `threshold`:         The misorientation threshold for being part of a family (deg)

    Returns the indices of the grain family
    """
    
    # Initialise
    lattice = crystallography.CubicLattice(1.0)
    crystal_direction = tensors.Vector(np.array([float(cd) for cd in crystal_direction])).normalize()
    sample_direction = tensors.Vector(np.array([float(sd) for sd in sample_direction])).normalize()
    rad_threshold = deg_to_rad(threshold)
    
    # Iterate through grains and add to family
    family_indices = []
    for i, orientation in enumerate(orientations):
        quat = rotations.Orientation(*orientation, convention="bunge")
        quat = quat.apply(sample_direction)
        dot_list = [crystal_direction.dot(sop.apply(quat)) for sop in lattice.symmetry.ops]
        misorientations = [np.arccos(np.abs(dot)) for dot in dot_list]
        misorientation = np.min(misorientations)
        if misorientation < rad_threshold:
            family_indices.append(i)
            
    # Return indices
    return family_indices

def transpose(list_of_lists:list) -> list:
    """
    Transposes a 2D list of lists
    
    Parameters:
    * `list_of_lists`: A list of lists (i.e., a 2D grid)
    
    Returns the transposed list of lists
    """
    transposed = np.array(list_of_lists).T.tolist()
    return transposed

def run_model(passive_eulers:list, weights:list) -> tuple:
    """
    Runs the crystal plasticity model
    
    Parameters:
    * `passive_euler`: List of euler-bunge angles (rads) in passive rotation
    * `weights`:       List of weights
    
    Returns the elastic model, single crystal model, polycrystal model, and driver results
    """
    
    # Get lattice
    lattice = crystallography.CubicLattice(1.0)
    lattice.add_slip_system([1,1,0], [1,1,1])

    # Get orientations
    active_eulers = [reorient_euler(euler) for euler in passive_eulers] # active
    orientations = [rotations.CrystalOrientation(*euler, angle_type="radians", convention="bunge") for euler in active_eulers]

    # Define elastic model
    # e_model = elasticity.IsotropicLinearElasticModel(211000, "youngs", 0.30, "poissons")
    # e_model = elasticity.CubicLinearElasticModel(284039, 121731, 162000, "components") # isotropic
    e_model = elasticity.CubicLinearElasticModel(205000, 138000, 126000, "components") # ondrej
    # e_model = elasticity.CubicLinearElasticModel(250000, 151000, 123000, "components") # k. mo
    # e_model = elasticity.CubicLinearElasticModel(170000, 121000, 75000, "components") # simmons
    # e_model = elasticity.CubicLinearElasticModel(210900, 140300, 122500, "components") # ledbetter
    
    # Define rest of model
    str_model  = slipharden.VoceSlipHardening(TAU_SAT, B, TAU_0)
    slip_model = sliprules.PowerLawSlipRule(str_model, GAMMA_0, N)
    i_model    = inelasticity.AsaroInelasticity(slip_model)
    k_model    = kinematics.StandardKinematicModel(e_model, i_model)
    sc_model   = singlecrystal.SingleCrystalModel(k_model, lattice, miter=16, max_divide=2, verbose=False)
    pc_model   = polycrystal.TaylorModel(sc_model, orientations, nthreads=NUM_THREADS, weights=weights)
    
    # Run driver and return everything
    results = drivers.uniaxial_test(pc_model, STRAIN_RATE, T=20, emax=MAX_STRAIN, nsteps=500, rtol=1e-6,
                                       atol=1e-10, miter=25, verbose=False, full_results=True)
    return e_model, sc_model, pc_model, results

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

# Define custom orientations
passive_eulers = deg_to_rad([[0,135,45], [68,38,72], [110,75,80], [0,0,0]])
weights = [1]*len(passive_eulers)

# # Define EBSD orientations
# grain_stats = np.loadtxt("../data/617_s3/grain_stats.csv", delimiter=",")[:-2]
# passive_eulers = [gs[:3] for gs in grain_stats]
# weights = [gs[3] for gs in grain_stats]

# Initialise
crystal_stress_list = [[] for _ in range(len(passive_eulers))]
crystal_strain_list = [[] for _ in range(len(passive_eulers))]

# Run the NEML model
e_model, sc_model, pc_model, results = run_model(passive_eulers, weights)
strain_list = np.array(results["strain"])[:,0]
stress_list = np.array(results["stress"])[:,0]
history     = np.array(results["history"])

# Plot stress-strain
plt.figure(figsize=(5,5))
plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
plt.xlabel("Strain (mm/mm)")
plt.ylabel("Stress (MPa)")
plt.plot(strain_list, stress_list)
plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white")
plt.savefig("plot_ss.png")
plt.cla()
plt.clf()

# Read stress tensors
nhist = sc_model.nstore
stress_xx = history[:, pc_model.n*nhist+0:pc_model.n*nhist+6*pc_model.n:6]
stress_yy = history[:, pc_model.n*nhist+1:pc_model.n*nhist+6*pc_model.n:6]
stress_zz = history[:, pc_model.n*nhist+2:pc_model.n*nhist+6*pc_model.n:6]
stress_xy = history[:, pc_model.n*nhist+3:pc_model.n*nhist+6*pc_model.n:6]
stress_yz = history[:, pc_model.n*nhist+4:pc_model.n*nhist+6*pc_model.n:6]
stress_xz = history[:, pc_model.n*nhist+5:pc_model.n*nhist+6*pc_model.n:6]

# from manual import orientation_to_c, multiply_tensors
# c_pi_list = [orientation_to_c(passive_euler) for passive_euler in passive_eulers]

# Calculate elastic strains
for i, state in enumerate(history):
    
    # Calculate NEML stress-strain
    orientation_list = pc_model.orientations(state)
    for j, orientation in enumerate(orientation_list):
        
        # Define stress tensor
        a = np.sqrt(2)
        stress_tensor = np.array([
            [stress_xx[i,j], stress_xy[i,j]/a, stress_xz[i,j]]/a,
            [stress_xy[i,j]/a, stress_yy[i,j], stress_yz[i,j]]/a,
            [stress_xz[i,j]/a, stress_yz[i,j]/a, stress_zz[i,j]]
        ])

        # strain_tensor = multiply_tensors(stress_tensor, c_pi_list[j])
        # crystal_stress_list[j].append(stress_tensor[0,0])
        # crystal_strain_list[j].append(strain_tensor[0,0]*1e6)

        # Define strain tensor
        active_orientation = orientation.inverse()
        # c_inverse = e_model.S_tensor(20, active_orientation)
        c_inverse = e_model.C_tensor(20, active_orientation).inverse()
        strain_product = c_inverse * tensors.Symmetric(stress_tensor)
        strain_tensor = strain_product.data
        
        # Store stress and strain values
        crystal_stress_list[j].append(stress_tensor[0,0])
        crystal_strain_list[j].append(strain_tensor[0]*1e6)

# Format plot
plt.figure(figsize=(5,5))
plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
plt.xlabel("Elastic Strain (μmm/mm)")
plt.ylabel("Applied Stress (MPa)")

# Iterate through the families
for family, colour in zip(FAMILY_LIST, COLOUR_LIST):
    
    # Identify family
    family_str = "".join([str(f) for f in family])
    family_indices = get_grain_family(passive_eulers, family, [1,0,0], 10)

    # Get strain and stress for family
    family_stress = [crystal_stress_list[i] for i in family_indices]
    family_stress = [np.average(fs) for fs in transpose(family_stress)]
    family_strain = [crystal_strain_list[i] for i in family_indices]
    family_strain = [np.average(fs) for fs in transpose(family_strain)]

    # Plot points
    plt.plot(family_strain, stress_list, color=colour, label=family_str)
    # plt.plot(family_strain, family_stress, color=colour, label=family_str)

# Save the plot
plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white")
plt.savefig("plot_es.png")

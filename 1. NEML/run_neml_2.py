"""
 Title:         Manual
 Description:   Manual calculation of elastic strain
 References:    [1] https://www.researchgate.net/publication/324088567_Computing_Euler_angles_with_Bunge_convention_from_rotation_matrix
 Author:        Janzen Choi

"""

# Libraries
import math, numpy as np
import matplotlib.pyplot as plt
from neml import elasticity, drivers
from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal
from neml.math import rotations, matrix, tensors
from pole_figure import IPF, get_lattice

# Define families
COLOUR_LIST = ["green", "black", "blue", "red"]
FAMILY_LIST = [[2,2,0], [1,1,1], [3,1,1], [2,0,0]]

# Constants
NUM_THREADS = 8
STRAIN_RATE = 1e-4
MAX_STRAIN  = 0.5

# Model parameters
TAU_0   = 112
GAMMA_0 = 3.3333e-5
N       = 15

def transpose(list_of_lists:list) -> list:
    """
    Transposes a 2D list of lists
    
    Parameters:
    * `list_of_lists`: A list of lists (i.e., a 2D grid)
    
    Returns the transposed list of lists
    """
    transposed = np.array(list_of_lists).T.tolist()
    return transposed

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

def run_model(passive_eulers:list, weights:list, sm_0:float, sm_1:float) -> tuple:
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

    # Define model
    e_model    = elasticity.CubicLinearElasticModel(250000, 151000, 123000, "components")
    # e_model    = elasticity.CubicLinearElasticModel(205000, 138000, 126000, "components") # ondrej
    
    sm_data    = [[sm_0]*lattice.ntotal for _ in range(lattice.ntotal)]
    for i in range(lattice.ntotal):
        sm_data[i][i] = sm_1
    sm_object  = matrix.SquareMatrix(lattice.ntotal, type="dense", data=np.array(sm_data).flatten())
    str_model  = slipharden.GeneralLinearHardening(sm_object, [TAU_0]*lattice.ntotal, absval=True)
    slip_model = sliprules.PowerLawSlipRule(str_model, GAMMA_0, N)
    
    i_model    = inelasticity.AsaroInelasticity(slip_model)
    k_model    = kinematics.StandardKinematicModel(e_model, i_model)
    sc_model   = singlecrystal.SingleCrystalModel(k_model, lattice, miter=16, max_divide=2, verbose=False)
    pc_model   = polycrystal.TaylorModel(sc_model, orientations, nthreads=NUM_THREADS, weights=weights)
    
    # Run driver and return everything
    results = drivers.uniaxial_test(pc_model, STRAIN_RATE, T=20, emax=MAX_STRAIN, nsteps=500, rtol=1e-6,
                                       atol=1e-10, miter=25, verbose=False, full_results=True)
    return e_model, sc_model, pc_model, results

# # Define custom orientations
# passive_eulers = deg_to_rad([[0,135,45], [68,38,72], [110,75,80], [0,0,0]])
# weights = [1]*len(passive_eulers)

# Define EBSD orientations
grain_stats    = np.loadtxt("grain_stats.csv", delimiter=",")[:5]
passive_eulers = [gs[:3] for gs in grain_stats]
weights        = [gs[3] for gs in grain_stats]

# Initialise IPF
ipf = IPF(get_lattice("fcc"))
direction = [1,0,0]
colour_list = ["red", "blue", "green", "orange", "purple"]

# Run the NEML model and plot results
for sm_0 in [1000]:
    for sm_1 in [2000]:
# for sm_0 in [5000*i for i in range(1,11)]:
#     for sm_1 in [5000*i for i in range(1,11)]:

        # Run the NEML model
        e_model, sc_model, pc_model, results = run_model(passive_eulers, weights, sm_0, sm_1)
        history = np.array(results["history"])

        strain_list = np.array(results["strain"])[:,0]
        stress_list = np.array(results["stress"])[:,0]

        # Plot stress-strain
        plt.figure(figsize=(5,5))
        plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
        plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":")
        plt.xlabel("Strain (mm/mm)")
        plt.ylabel("Stress (MPa)")
        plt.plot(strain_list, stress_list)
        plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white")
        plt.savefig("plot_ss.png")
        plt.cla()
        plt.clf()
        exit()

        # # Calculate reorientation trajectories
        # orientation_history = [[list(orientation.to_euler(angle_type="radians", convention="bunge"))
        #                         for orientation in pc_model.orientations(state)] for state in history]
        # trajectories = [[] for _ in range(len(passive_eulers))]
        # for state in history:
        #     for i, orientation in enumerate(pc_model.orientations(state)):
        #         euler = list(orientation.to_euler(angle_type="radians", convention="bunge"))
        #         trajectories[i].append(euler)
                
        # # Plot the IPF
        # for colour, trajectory in zip(colour_list, trajectories):
        #     ipf.plot_ipf_trajectory([trajectory], direction, "plot", {"color": colour, "linewidth": 1, "zorder": 3})
        #     ipf.plot_ipf_trajectory([trajectory], direction, "arrow", {"color": colour, "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
        #     ipf.plot_ipf_trajectory([[trajectory[0]]], direction, "scatter", {"color": colour, "s": 6**2, "zorder": 3})
        # plt.savefig("plot_rt.png")

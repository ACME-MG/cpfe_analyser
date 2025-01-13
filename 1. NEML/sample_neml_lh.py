"""
 Title:         Manual
 Description:   Manual calculation of elastic strain
 References:    [1] https://www.researchgate.net/publication/324088567_Computing_Euler_angles_with_Bunge_convention_from_rotation_matrix
                [2] https://pdf.sciencedirectassets.com/271567/1-s2.0-S0749641913X00026/1-s2.0-S0749641912001635/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEPr%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIBd5u6Hel8h21ZBNTE3ssZC3U4k2dNfF%2BlrT96Yfz3J7AiEAryt%2FKXUre7Vr3KC06OI7xTRtS%2BITfcm%2BZ1iYkwI%2BjkIqsgUIUhAFGgwwNTkwMDM1NDY4NjUiDPxCAb2vsQKyAgBjmyqPBX%2FnPzY3Nzj8Q4H0bVggHkPZNoojAyM22e7zyhWNk29EHkqeYbEH6cOGXfTBzDR4Lm962gI51lwuQXBB7vC7NIPBs9FiF3Z6O2HdRybHpJ80dFEgYKySeHSCg5YccjEXbFNcnXDdeYHFpnFQdl7n%2B3uJ9of1tz1ifOYHzBgqBTWGjUyJUdcR1DDxGguzxEMGf4JZOCM26bK0oQdyD0iXOfW5XfRg65aa34JbTiuEj4OZa%2FHzsna482j9Rgo1%2FYUXPp9wufWlEuqWjFxEeog%2BBIzX9unVQyuz%2BfZsjHjZ3VhLzBubyg%2BbOrKsudYNSOeEXS2z3KnCxxdFtX3mEmLbSg%2FynKm7SooLMMSk248m1D12x3iH3tPjYnW6vDMc%2FvY2V1srWm1MNlmXJSJ4SjkT31fHZGMHqEErPTmL3VZOI%2FxsE02Bh%2FNlIy7aFRKTbrKBYMB%2BLegbhBAqj6eJYzHbavs85aQ2Glyf7mQS%2B2W49XOzPRpaOmeMjkbCYSJOUnbc%2Fk5RLrVYPhQ8lbjsd7VhmPVxdBAWELPS%2BokMXYfIsEsD5abYFxOhg1X9XrB6Pc5Mj8EjZ5zB%2BqSXP23dxc0CQTi8CjqzdJh1uN37VscGwnh07h6fIZxKgNMSq%2B1kksHJnaHDGb7zqfHwOXL0Bl4vgqvdKrfvUnVXm95LgioU8JYi6aEGrphdJm8RGDvvqvTsYqwfSb4WOXB14HizNWPpYSMqZHzxR9McAMfTLce5eKwPHbjvW5449qHQnvDZIYOzgyYi13E8jJh7nqUqcHyN2w%2BsmOOgTkQoVs%2Fy0MNWr0xszY367EE8IccUVE6l3fCQsyJ%2BfKY1yPhHdyLRBk1aHAP5FsKieNX8J9SN%2BYtIKTgwtrGXuAY6sQEyENrWgFSOjkPwzuNgiFTFk%2B563V402Ya%2F%2FLHQVztH0QrFlV4u0TL07GSCtInZeT%2FlcMG5I%2FVtMFSXgDi1mq%2BnqsUjZ%2FAyb9Fi4Ei4yQINtAiWwhelY%2B1N65jbCGoT9Mp%2BatxfocRTeZeOGnhGwaszVbmEcOnXjWSO%2BEqqSZTMVrCOPd0lBkE1vp8COKM4aGhTyTGrEFPZP0iLpM870CIt%2B8xzOFaRHGujaFqNaf587Yc%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20241009T015940Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYUQF4X53J%2F20241009%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=ea74c3b5f90a6ec9fd47c08b35dd2234fd4694044b5a3cf8952beb40fbf9224c&hash=51d2311a9ba522febc66efa72c488a6cc5df37bfef658c6e52017d7bcd6add84&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0749641912001635&tid=spdf-8e330e4f-a442-464f-b587-27a11e861216&sid=57902e0d439955488648e384d345c8bf9647gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=07135b0301030057005e&rr=8cfac5c82c86a817&cc=au
 Author:        Janzen Choi

"""

# Libraries
import math, numpy as np
import matplotlib.pyplot as plt
import sys; sys.path += [".."]
import pyDOE2
from neml import elasticity, drivers
from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal
from neml.math import rotations, matrix, tensors
from __common__.pole_figure import IPF, get_lattice

# Define families
COLOUR_LIST = ["green", "black", "blue", "red"]
FAMILY_LIST = [[2,2,0], [1,1,1], [3,1,1], [2,0,0]]

# Constants
NUM_THREADS = 8
STRAIN_RATE = 1e-4
MAX_STRAIN  = 0.3

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

def linear_scale(value:float, in_l_bound:float, in_u_bound:float, out_l_bound:float, out_u_bound:float) -> float:
    """
    Linearly scales a value

    Parameters:
    * `value`:       The value to be linearly scaled
    * `in_l_bound`:  The lower bound of the input values
    * `in_u_bound`:  The upper bound of the input values
    * `out_l_bound`: The lower bound of the linearly scaled values
    * `out_u_bound`: The upper bound of the linearly scaled values
    """
    in_range = in_u_bound - in_l_bound
    out_range = out_u_bound - out_l_bound
    scaled_value = (value-in_l_bound)*out_range/in_range + out_l_bound
    return scaled_value

def get_lhs(param_bounds:dict, num_points:int) -> list:
    """
    Generates points using latin hypercube sampling
    
    Parameters:
    * `param_bounds`: Dictionary of parameter bounds;
                      (i.e., `{"name": (l_bound, u_bound)}`)
    * `num_points`:   The number of points to sample
    
    Returns the list of dictionaries of parameter combinations
    """
    
    # Get unscaled LHS points
    num_params = len(param_bounds.keys())
    combinations = list(pyDOE2.lhs(num_params, samples=num_points))
    
    # Create linear scales for each parameter
    param_scales = {}
    for param in param_bounds.keys():
        param_scales[param] = lambda x : linear_scale(x, 0, 1, param_bounds[param][0], param_bounds[param][1])

    # Linearly scale the unscaled combinations
    scaled_dict_list = []
    for combination in combinations:
        scaled_dict = {}
        for i, param in enumerate(param_scales.keys()):
            scaled_param = param_scales[param](combination[i])
            scaled_dict[param] = scaled_param
        scaled_dict_list.append(scaled_dict)

    # Return scaled LHS points
    return scaled_dict_list

def get_big_lh_matrix(params:tuple) -> list:
    """
    Gets the 12x12 latent hardening matrix
    
    Parameters:
    * `params`: Elements of the matrix
    
    Returns the matrix in the form of a list of lists
    """
    h_0, h_1, h_2, h_3, h_4, h_5 = params
    lh_matrix = [
        [h_0, h_1, h_1, h_3, h_4, h_4, h_2, h_4, h_5, h_2, h_5, h_4],
        [h_1, h_0, h_1, h_4, h_2, h_5, h_4, h_3, h_4, h_5, h_2, h_4],
        [h_1, h_1, h_0, h_4, h_5, h_2, h_5, h_4, h_2, h_4, h_4, h_3],
        [h_3, h_4, h_4, h_0, h_1, h_1, h_2, h_5, h_4, h_2, h_4, h_5],
        [h_4, h_2, h_5, h_1, h_0, h_1, h_5, h_2, h_4, h_4, h_3, h_4],
        [h_4, h_5, h_2, h_1, h_1, h_0, h_4, h_4, h_3, h_5, h_4, h_2],
        [h_2, h_4, h_5, h_2, h_5, h_4, h_0, h_1, h_1, h_3, h_4, h_4],
        [h_4, h_3, h_4, h_5, h_2, h_4, h_1, h_0, h_1, h_4, h_2, h_5],
        [h_5, h_4, h_2, h_4, h_4, h_3, h_1, h_1, h_0, h_4, h_5, h_2],
        [h_2, h_5, h_4, h_2, h_4, h_5, h_3, h_4, h_4, h_0, h_1, h_1],
        [h_5, h_2, h_4, h_4, h_3, h_4, h_4, h_2, h_5, h_1, h_0, h_1],
        [h_4, h_4, h_3, h_5, h_4, h_2, h_4, h_5, h_2, h_1, h_1, h_0]
    ]
    return lh_matrix

def get_lh_matrix(params:tuple) -> list:
    """
    Gets the 12x12 latent hardening matrix
    
    Parameters:
    * `params`: Elements of the matrix
    
    Returns the matrix in the form of a list of lists
    """
    h_0, h_1 = params
    lh_matrix = [[h_0]*12 for _ in range(12)]
    for i in range(12):
        lh_matrix[i][i] = h_1
    return lh_matrix

def csv_to_dict(csv_path:str, delimeter:str=",") -> dict:
    """
    Converts a CSV file into a dictionary
    
    Parameters:
    * `csv_path`:  The path to the CSV file
    * `delimeter`: The separating character
    
    Returns the dictionary
    """

    # Read all data from CSV (assume that file is not too big)
    csv_fh = open(csv_path, "r", encoding="utf-8-sig")
    csv_lines = csv_fh.readlines()
    csv_fh.close()

    # Initialisation for conversion
    csv_dict = {}
    headers = csv_lines[0].replace("\n", "").split(delimeter)
    csv_lines = csv_lines[1:]
    for header in headers:
        csv_dict[header] = []

    # Start conversion to dict
    for csv_line in csv_lines:
        csv_line_list = csv_line.replace("\n", "").split(delimeter)
        for i in range(len(headers)):
            value = csv_line_list[i]
            if value == "":
                continue
            try:
                value = float(value)
            except:
                pass
            csv_dict[headers[i]].append(value)
    
    # Convert single item lists to items and things multi-item lists
    for header in headers:
        if len(csv_dict[header]) == 1:
            csv_dict[header] = csv_dict[header][0]
    
    # Return
    return csv_dict

def run_model(tau_0:float, gamma_0:float, n:float, lh_0:float, lh_1:float,
            #   lh_2:float, lh_3:float, lh_4:float, lh_5:float,
              passive_eulers:list, weights:list) -> tuple:
    """
    Runs the crystal plasticity model
    
    Parameters:
    * `passive_euler`: List of euler-bunge angles (rads) in passive rotation
    * `weights`:       List of weights
    
    Returns the elastic model, single crystal model, polycrystal model, and driver results
    """
    
    # Change LH parameters
    # lh_1 *= lh_0
    # lh_2 *= lh_0
    # lh_3 *= lh_0
    # lh_4 *= lh_0
    # lh_5 *= lh_0

    # Get lattice
    lattice = crystallography.CubicLattice(1.0)
    lattice.add_slip_system([1,1,0], [1,1,1])

    # Get orientations
    active_eulers = [reorient_euler(euler) for euler in passive_eulers] # active
    orientations = [rotations.CrystalOrientation(*euler, angle_type="radians", convention="bunge") for euler in active_eulers]

    # Define model
    e_model    = elasticity.CubicLinearElasticModel(250000, 151000, 123000, "components")
    # lh_matrix  = get_big_lh_matrix((lh_0, lh_1, lh_2, lh_3, lh_4, lh_5))
    lh_matrix  = get_lh_matrix((lh_0, lh_1))
    sm_object  = matrix.SquareMatrix(lattice.ntotal, type="dense", data=np.array(lh_matrix).flatten())
    str_model  = slipharden.GeneralLinearHardening(sm_object, [tau_0]*lattice.ntotal, absval=True)
    slip_model = sliprules.PowerLawSlipRule(str_model, gamma_0, n)
    i_model    = inelasticity.AsaroInelasticity(slip_model)
    k_model    = kinematics.StandardKinematicModel(e_model, i_model)
    sc_model   = singlecrystal.SingleCrystalModel(k_model, lattice, miter=16, max_divide=2, verbose=False)
    pc_model   = polycrystal.TaylorModel(sc_model, orientations, nthreads=NUM_THREADS, weights=weights)
    
    # Run driver and return everything
    results = drivers.uniaxial_test(pc_model, STRAIN_RATE, T=20, emax=MAX_STRAIN, nsteps=100, rtol=1e-6,
                                       atol=1e-10, miter=25, verbose=False, full_results=True)
    return e_model, sc_model, pc_model, results

# Read experimental data
exp_data = csv_to_dict("data/617_s3_exp.csv")
exp_grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in exp_data.keys() if "_phi_1" in key][:5]
exp_trajectories = [transpose([exp_data[f"g{grain_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]])
                    for grain_id in exp_grain_ids]

# Define EBSD orientations
grain_stats    = np.loadtxt("data/grain_stats.csv", delimiter=",")[:10]
passive_eulers = [gs[:3] for gs in grain_stats]
weights        = [gs[3] for gs in grain_stats]

# Initialise stress-strain plot
plt.figure(1, figsize=(5,5))
plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
plt.xlabel("Strain (mm/mm)")
plt.ylabel("Stress (MPa)")
plt.scatter(exp_data["strain"], exp_data["stress"], color="silver", s=8**2)

# Initialise IPF
plt.figure(2)
ipf = IPF(get_lattice("fcc"))
direction = [1,0,0]
colour_list = ["red", "blue", "green", "orange", "purple"]

# Plot trajectories
ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": "silver", "linewidth": 2})
ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": "silver", "head_width": 0.01, "head_length": 0.015})
ipf.plot_ipf_trajectory([[et[0]] for et in exp_trajectories], direction, "scatter", {"color": "silver", "s": 8**2})
for exp_trajectory, grain_id in zip(exp_trajectories, exp_grain_ids):
    ipf.plot_ipf_trajectory([[exp_trajectory[0]]], direction, "text", {"color": "black", "fontsize": 8, "s": grain_id})

# Get CP parameter combinations
bounds_dict = {
    "cp_lh_0":    (0, 400), # main diagonal
    "cp_lh_1":    (0, 400),
    # "cp_lh_2":    (0, 800),
    # "cp_lh_3":    (0, 800),
    # "cp_lh_4":    (0, 800),
    # "cp_lh_5":    (0, 800),
    "cp_tau_0":   (0, 200), # yield/3
    "cp_n":       (1, 16),
    "cp_gamma_0": (1e-4/3, 1e-4/3),
    # "cp_tau_0":   (0, 200), # yield/3
    # "cp_n":       (1, 16),
    # "cp_gamma_0": (1e-4/3, 1e-4/3)
}
params_list = get_lhs(bounds_dict, 32)
params_list = [dict(zip(bounds_dict.keys(), [398.21, 6.856, 69.407, 15.658] + [1e-4/3]))]

# Iterate through parameters
for count, params in enumerate(params_list):

    # Print progress
    print(f"===============[{count}]===============")
    for param_name in params.keys():
        print(f"{param_name}:\t{params[param_name]}")
    print(f"===============[{count}]===============")
    
    # Run the NEML model
    try:
        e_model, sc_model, pc_model, results = run_model(
            tau_0          = params["cp_tau_0"],
            gamma_0        = params["cp_gamma_0"],
            n              = params["cp_n"],
            lh_0           = params["cp_lh_0"],
            lh_1           = params["cp_lh_1"],
            # lh_2           = params["cp_lh_2"],
            # lh_3           = params["cp_lh_3"],
            # lh_4           = params["cp_lh_4"],
            # lh_5           = params["cp_lh_5"],
            passive_eulers = passive_eulers,
            weights        = weights
        )
    except:
        print("\n[Failed]\n")
        continue
    
    # Extract data
    history = np.array(results["history"])
    strain_list = np.array(results["strain"])[:,0]
    stress_list = np.array(results["stress"])[:,0]

    # Plot stress-strain
    plt.figure(1)
    plt.plot(strain_list, stress_list)
    plt.text(strain_list[-1], stress_list[-1], f"{count}")
    plt.savefig("results/plot_ss.png")

    # Calculate reorientation trajectories
    orientation_history = [[list(orientation.to_euler(angle_type="radians", convention="bunge"))
                            for orientation in pc_model.orientations(state)] for state in history]
    trajectories = [[] for _ in range(len(passive_eulers))]
    for state in history:
        for i, orientation in enumerate(pc_model.orientations(state)):
            euler = list(orientation.to_euler(angle_type="radians", convention="bunge"))
            trajectories[i].append(euler)
            
    # Plot the IPF
    plt.figure(2)
    for colour, trajectory in zip(colour_list, trajectories):
        ipf.plot_ipf_trajectory([trajectory], direction, "plot", {"color": colour, "linewidth": 1, "zorder": 3})
        ipf.plot_ipf_trajectory([trajectory], direction, "arrow", {"color": colour, "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
        ipf.plot_ipf_trajectory([[trajectory[0]]], direction, "scatter", {"color": colour, "s": 6**2, "zorder": 3})
    plt.savefig("results/plot_rt.png")

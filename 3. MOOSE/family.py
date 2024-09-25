"""
 Title:         Family
 Description:   For testing grain family grouping
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += ["/home/janzen/code/moose_sim"]
from moose_sim.analyse.plotter import save_plot
from moose_sim.analyse.pole_figure import IPF, get_lattice
from moose_sim.helper.general import transpose, remove_consecutive_duplicates
from moose_sim.helper.io import csv_to_dict
from moose_sim.maths.familiser import get_grain_family

def get_trajectories(data_dict:dict) -> list:
    """
    Gets the reorientation trajectories

    Parameters:
    * `data_dict`: The data dictionary
    
    Return trajectories as a list of lists of euler angles
    """
    grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in data_dict.keys() if "_phi_1" in key]
    trajectories = []
    for grain_id in grain_ids:
        trajectory = [data_dict[f"g{grain_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]
        trajectory = transpose(trajectory)
        trajectory = remove_consecutive_duplicates(trajectory)
        trajectories.append(trajectory)
    return trajectories

# Initialise IPF plotter
ipf = IPF(get_lattice("fcc"))#, sample_symmetry=[1])
sample_direction = [1,0,0]
crystal_direction = [2,2,0]

# Get experimental orientations
exp_dict = csv_to_dict("../data/617_s3/617_s3_exp.csv")
exp_orientations = [trajectory[0] for trajectory in get_trajectories(exp_dict)]
ipf.plot_ipf_trajectory([[eo] for eo in exp_orientations], sample_direction, "scatter", {"color": "darkgray", "s": 8**2})

# Plot grain family on IPF
family_indices = get_grain_family(exp_orientations, crystal_direction, sample_direction, 10)
exp_orientations = [exp_orientations[i] for i in family_indices]
ipf.plot_ipf_trajectory([[eo] for eo in exp_orientations], sample_direction, "scatter", {"color": "red", "s": 5**2})
save_plot("results/family.png")

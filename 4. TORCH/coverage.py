"""
 Title:         Coverage
 Description:   Plots multiple simulation results
 Author:        Janzen Choi

"""

# Libraries
import os
import sys; sys.path += [".."]
from __common__.io import csv_to_dict
from __common__.general import transpose
from __common__.pole_figure import get_lattice, IPF
from __common__.plotter import save_plot, Plotter

# Constants
EXP_PATH     = "data/617_s3_exp.csv"
RESULTS_DIR  = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim"
SAMPLED_PATH = f"{RESULTS_DIR}/2025-01-07 (617_s3_40um_lh_sm32)"
# STRAIN_FIELD = "average_grain_strain"
# STRESS_FIELD = "average_grain_stress"
STRAIN_FIELD = "average_strain"
STRESS_FIELD = "average_stress"
RT_COLOUR    = "green"

# Main function
def main():

    # Get experimental data
    exp_dict = csv_to_dict(EXP_PATH)

    # Get simulated data
    dir_path_list = [f"{SAMPLED_PATH}/{dir_path}" for dir_path in os.listdir(SAMPLED_PATH)
                    if os.path.exists(f"{SAMPLED_PATH}/{dir_path}/summary.csv")]
    sim_path_list = [f"{dir_path}/summary.csv" for dir_path in dir_path_list]
    sim_dict_list = [csv_to_dict(sim_path) for sim_path in sim_path_list]

    # Plot stress-strain coverage
    plot_ss(exp_dict, sim_dict_list, "results/plot_cvg_ss.png")

    # Plot reorientation trajectory coverage
    all_grain_ids = sorted([int(key.replace("_phi_1","").replace("g","")) for key in sim_dict_list[0].keys() if "_phi_1" in key])
    all_grain_ids = [59, 63, 86, 237, 303]
    num_splits = 5

    grain_ids_list = [all_grain_ids[i:i + num_splits] for i in range(0, len(all_grain_ids), num_splits)]
    for i, grain_ids in enumerate(grain_ids_list):
        plot_rt(exp_dict, sim_dict_list, grain_ids, f"results/plot_cvg_rt_{i+1}")

def plot_ss(exp_dict:dict, sim_dict_list:list, path:str=""):
    """
    Plots the stress-strain response of the sampled simulations
    
    Parameters:
    * `exp_dict`:      Dictionary of experimental data
    * `sim_dict_list`: List of dictionaries of simulation results
    * `path`:          The path to save the plot
    """
    plotter = Plotter("strain", "stress", "mm/mm", "MPa")
    plotter.prep_plot()
    plotter.scat_plot(exp_dict, "silver", "Experimental")
    for sim_dict in sim_dict_list:
        sim_dict["strain"] = sim_dict[STRAIN_FIELD]
        sim_dict["stress"] = sim_dict[STRESS_FIELD]
        plotter.line_plot(sim_dict, "green", "Calibration")
    save_plot(path)

def plot_rt(exp_dict:dict, sim_dict_list:list, grain_ids:list, path:str=""):
    """
    Plots the reorientation trajectories response of the sampled simulations
    
    Parameters:
    * `exp_dict`:      Dictionary of experimental data
    * `sim_dict_list`: List of dictionaries of simulation results
    * `grain_ids`: List of grain IDs
    * `path`:      The path to save the plot
    """

    # Initialise IPF
    ipf = IPF(get_lattice("fcc"))
    direction = [1,0,0]
    get_trajectories = lambda dict : [transpose([dict[f"g{grain_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]) for grain_id in grain_ids]
    
    # Plot experimental reorientation trajectories
    exp_trajectories = get_trajectories(exp_dict)
    ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": "silver", "linewidth": 2})
    ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": "silver", "head_width": 0.01, "head_length": 0.015})
    ipf.plot_ipf_trajectory([[et[0]] for et in exp_trajectories], direction, "scatter", {"color": "silver", "s": 8**2})

    # Plot simulated reorientation trajectories
    for sim_dict in sim_dict_list:
        cal_trajectories = get_trajectories(sim_dict)
        ipf.plot_ipf_trajectory(cal_trajectories, direction, "plot", {"color": RT_COLOUR, "linewidth": 1, "zorder": 3})
        ipf.plot_ipf_trajectory(cal_trajectories, direction, "arrow", {"color": RT_COLOUR, "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
        ipf.plot_ipf_trajectory([[ct[0]] for ct in cal_trajectories], direction, "scatter", {"color": RT_COLOUR, "s": 6**2, "zorder": 3})
    
    # Plot grain ID
    for exp_trajectory, grain_id in zip(exp_trajectories, grain_ids):
        ipf.plot_ipf_trajectory([[exp_trajectory[0]]], direction, "text", {"color": "blue", "fontsize": 8, "s": grain_id, "zorder": 3})

    # Save IPF
    save_plot(path)

# Calls the main function
if __name__ == "__main__":
    main()

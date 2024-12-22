"""
 Title:         Plot Optimised
 Description:   Plots the response of the optimised simulation
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
from __common__.io import csv_to_dict
from __common__.general import transpose
from __common__.pole_figure import get_lattice, IPF
from __common__.plotter import define_legend, save_plot, Plotter

# Constants
EXP_PATH = "data/617_s3_exp.csv"
# SIM_FILE = "2024-11-05 (617_s3_40um_lh2_opt)"
SIM_FILE = "2024-12-22 (max_0p3_i12)/241222090646_i12_simulate"
SIM_PATH = f"/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo/{SIM_FILE}/summary.csv"
CAL_GRAIN_IDS = [59, 63, 86, 237, 303]
VAL_GRAIN_IDS = [44, 53, 60, 78, 190]
# STRAIN_FIELD = "average_grain_strain"
# STRESS_FIELD = "average_grain_stress"
STRAIN_FIELD = "average_strain"
STRESS_FIELD = "average_stress"

# Main function
def main():

    # Get all results
    res_dict = csv_to_dict(SIM_PATH)
    exp_dict = csv_to_dict(EXP_PATH)

    # Plot reorientation trajectories
    plot_trajectories(exp_dict, res_dict, CAL_GRAIN_IDS, "green", "Calibration", "results/plot_opt_cal_rt.png")
    plot_trajectories(exp_dict, res_dict, VAL_GRAIN_IDS, "red",   "Validation",  "results/plot_opt_val_rt.png")

    # Plot stress-strain curve
    res_dict["strain"] = res_dict[STRAIN_FIELD]
    res_dict["stress"] = res_dict[STRESS_FIELD]
    plotter = Plotter("strain", "stress", "mm/mm", "MPa")
    plotter.prep_plot()
    plotter.scat_plot(exp_dict, "silver", "Experimental")
    plotter.line_plot(res_dict, "green", "Calibration")
    plotter.set_legend()
    save_plot("results/plot_opt_ss.png")

def plot_trajectories(exp_dict:dict, sim_dict:dict, grain_ids:list, sim_colour:str,
                      sim_label:str, path:str) -> None:
    """
    Plots the experimental and simulated reorientation trajectories

    Parameters:
    * `exp_dict`:   Dictionary of experimental data
    * `sim_dict`:   Dictionary of simulated data
    * `grain_ids`:  List of grain IDs
    * `sim_colour`: Colour to plot the simulated trajectories
    * `sim_label`:  Label for the simulated data
    * `path`:       Path to save the plot
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
    for exp_trajectory, grain_id in zip(exp_trajectories, grain_ids):
        ipf.plot_ipf_trajectory([[exp_trajectory[0]]], direction, "text", {"color": "black", "fontsize": 8, "s": grain_id})

    # Plot calibration reorientation trajectories
    sim_trajectories = get_trajectories(sim_dict)
    ipf.plot_ipf_trajectory(sim_trajectories, direction, "plot", {"color": sim_colour, "linewidth": 1, "zorder": 3})
    ipf.plot_ipf_trajectory(sim_trajectories, direction, "arrow", {"color": sim_colour, "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
    ipf.plot_ipf_trajectory([[st[0]] for st in sim_trajectories], direction, "scatter", {"color": sim_colour, "s": 6**2, "zorder": 3})

    # Save IPF
    define_legend(["silver", sim_colour], ["Experimental", sim_label], ["scatter", "line"])
    save_plot(path)

# Calls the main function
if __name__ == "__main__":
    main()

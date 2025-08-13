"""
 Title:         Plot Optimised Simulation
 Description:   Plots the response of the optimised simulation
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import matplotlib.pyplot as plt
from __common__.io import csv_to_dict
from __common__.general import transpose
from __common__.pole_figure import get_lattice, IPF
from __common__.plotter import save_plot, Plotter

# Paths
EXP_PATH = "data/617_s3_exp.csv"
SIM_FILE = "2025-01-18 (lh6_sm72_i20)/250117013234_i11_simulate"
SIM_PATH = f"/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo/{SIM_FILE}/summary.csv"
# SIM_FILE = "2025-01-01 (617_s3_10um_vh)"
# SIM_PATH = f"/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/moose_sim/{SIM_FILE}/summary.csv"

# Simulation information
CAL_GRAIN_IDS = [14, 72, 95, 101, 207, 240, 262, 287]
VAL_GRAIN_IDS = [39, 50, 138, 164, 185, 223, 243, 238]
STRAIN_FIELD = "average_strain"
STRESS_FIELD = "average_stress"
RES_DATA_MAP = "data/res_grain_map.csv"

# Initial EBSD maps
EXP_EBSD_ID = "ebsd_4"
# SIM_EBSD_ID = "ebsd_2"
SIM_EBSD_ID = "ebsd_4"

# Colours
EXP_COLOUR = "silver"
CAL_COLOUR = "tab:green"
VAL_COLOUR = "tab:red"

# Toggle parameters
SHOW_GRAIN_ID = False

# Main function
def main():

    # Get all results
    res_dict = csv_to_dict(SIM_PATH)
    exp_dict = csv_to_dict(EXP_PATH)

    # Plot reorientation trajectories
    plot_trajectories(exp_dict, res_dict, CAL_GRAIN_IDS, CAL_COLOUR, "Calibration", "results/plot_opt_cal_rt.png")
    plot_trajectories(exp_dict, res_dict, VAL_GRAIN_IDS, VAL_COLOUR, "Validation",  "results/plot_opt_val_rt.png")

    # Plot stress-strain curve
    res_dict["strain"] = res_dict[STRAIN_FIELD]
    res_dict["stress"] = res_dict[STRESS_FIELD]
    plotter = Plotter("strain", "stress", "mm/mm", "MPa")
    plotter.prep_plot(size=14)
    plotter.set_limits((0,0.5), (0,1400))
    plt.scatter(exp_dict["strain"], exp_dict["stress"], color=EXP_COLOUR, s=8**2)
    plt.plot(res_dict["strain"], res_dict["stress"], color=CAL_COLOUR, linewidth=3)
    handles = [
        plt.scatter([], [], color=EXP_COLOUR, label="Experiment", s=8**2),
        plt.plot([], [],    color=CAL_COLOUR, label="Calibration", linewidth=3)[0],
    ]
    legend = plt.legend(handles=handles, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    plt.gca().add_artist(legend)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
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
    get_trajectories = lambda dict, g_ids : [transpose([dict[f"g{g_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]) for g_id in g_ids]

    # Plot experimental reorientation trajectories
    exp_trajectories = get_trajectories(exp_dict, grain_ids)
    ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": EXP_COLOUR, "linewidth": 3})
    ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": EXP_COLOUR, "head_width": 0.01, "head_length": 0.015})
    ipf.plot_ipf_trajectory([[et[0]] for et in exp_trajectories], direction, "scatter", {"color": EXP_COLOUR, "s": 8**2})

    # Plot simulated reorientation trajectories
    sim_grain_ids = [get_sim_grain_id(grain_id) for grain_id in grain_ids]
    sim_trajectories = get_trajectories(sim_dict, sim_grain_ids)
    if 44 in sim_grain_ids:
        sim_trajectories[sim_grain_ids.index(44)] = sim_trajectories[sim_grain_ids.index(44)][:-9] 
    ipf.plot_ipf_trajectory(sim_trajectories, direction, "plot", {"color": sim_colour, "linewidth": 2, "zorder": 3})
    ipf.plot_ipf_trajectory(sim_trajectories, direction, "arrow", {"color": sim_colour, "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
    ipf.plot_ipf_trajectory([[st[0]] for st in sim_trajectories], direction, "scatter", {"color": sim_colour, "s": 6**2, "zorder": 3})

    # Plot grain IDs
    if SHOW_GRAIN_ID:
        for exp_trajectory, grain_id in zip(exp_trajectories, grain_ids):
            ipf.plot_ipf_trajectory([[exp_trajectory[0]]], direction, "text", {"color": "blue", "fontsize": 8, "s": grain_id, "zorder": 3})
    
    # Format and save IPF plot
    handles = [
        plt.plot([], [], color=EXP_COLOUR, label="Experiment", linewidth=3)[0],
        plt.plot([], [], color=sim_colour, label=sim_label,    linewidth=2)[0]
    ]
    legend = plt.legend(handles=handles, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    plt.gca().add_artist(legend)
    save_plot(path)

def get_sim_grain_id(exp_grain_id:int) -> int:
    """
    Maps the experimental to simulated grain ID

    Parameters:
    * `exp_grain_id`: The grain ID from the experimental data

    Returns the corresponding grain ID in the simulated data
    """
    grain_map = csv_to_dict(RES_DATA_MAP)
    exp_ebsd_ids = grain_map[EXP_EBSD_ID]
    sim_ebsd_ids = grain_map[SIM_EBSD_ID]
    return int(sim_ebsd_ids[exp_ebsd_ids.index(exp_grain_id)])

# Calls the main function
if __name__ == "__main__":
    main()

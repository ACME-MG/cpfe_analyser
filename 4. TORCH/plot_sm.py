"""
 Title:         Plot Surrogate
 Description:   Plots the response of the surrogate model
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += ["..", "/home/janzen/code/mms"]
from __common__.io import csv_to_dict
from __common__.general import transpose
from __common__.pole_figure import get_lattice, IPF
from __common__.plotter import define_legend, save_plot, Plotter
from __common__.surrogate import Model

# Paths
EXP_PATH = "data/617_s3_exp.csv"
# DIRECTORY = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo/2024-12-12 test_run/241213172919_i1_s1_sm"
DIRECTORY = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/mms/2024-12-13 (617_s3_40um_lh2_s32)"
SUR_PATH = f"{DIRECTORY}/sm.pt"
MAP_PATH = f"{DIRECTORY}/map.csv"

# Constants
MAX_STRAIN = 0.1
PARAM_LIST = [10.607, 120.94, 396.76, 1.0305, 7.26E-05]
CAL_GRAIN_IDS = [59, 63, 86, 237, 303]
VAL_GRAIN_IDS = [44, 53, 60, 78, 190]

# Main function
def main():

    # Get all results
    model = Model(SUR_PATH, MAP_PATH, EXP_PATH, MAX_STRAIN)
    res_dict = model.get_response(PARAM_LIST)
    exp_dict = csv_to_dict(EXP_PATH)
    get_trajectories = lambda dict, grain_ids : [transpose([dict[f"g{grain_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]) for grain_id in grain_ids]

    # Initialise IPF
    ipf = IPF(get_lattice("fcc"))
    direction = [1,0,0]

    # Plot experimental reorientation trajectories
    exp_trajectories = get_trajectories(exp_dict, CAL_GRAIN_IDS+VAL_GRAIN_IDS)
    ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": "silver", "linewidth": 2})
    ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": "silver", "head_width": 0.01, "head_length": 0.015})
    ipf.plot_ipf_trajectory([[et[0]] for et in exp_trajectories], direction, "scatter", {"color": "silver", "s": 8**2})
    for exp_trajectory, grain_id in zip(exp_trajectories, CAL_GRAIN_IDS+VAL_GRAIN_IDS):
        ipf.plot_ipf_trajectory([[exp_trajectory[0]]], direction, "text", {"color": "black", "fontsize": 8, "s": grain_id})

    # Plot calibration reorientation trajectories
    cal_trajectories = get_trajectories(res_dict, CAL_GRAIN_IDS)
    ipf.plot_ipf_trajectory(cal_trajectories, direction, "plot", {"color": "green", "linewidth": 1, "zorder": 3})
    ipf.plot_ipf_trajectory(cal_trajectories, direction, "arrow", {"color": "green", "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
    ipf.plot_ipf_trajectory([[ct[0]] for ct in cal_trajectories], direction, "scatter", {"color": "green", "s": 6**2, "zorder": 3})

    # Plot calibration reorientation trajectories
    val_trajectories = get_trajectories(res_dict, VAL_GRAIN_IDS)
    ipf.plot_ipf_trajectory(val_trajectories, direction, "plot", {"color": "red", "linewidth": 1, "zorder": 3})
    ipf.plot_ipf_trajectory(val_trajectories, direction, "arrow", {"color": "red", "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
    ipf.plot_ipf_trajectory([[vt[0]] for vt in val_trajectories], direction, "scatter", {"color": "red", "s": 6**2, "zorder": 3})

    # Save IPF
    define_legend(["silver", "green", "red"], ["Experimental", "Calibration", "Validation"], ["scatter", "line", "line"])
    save_plot("results/plot_sm_rt.png")

    # Plot stress-strain curve
    if "stress" in res_dict.keys():
        plotter = Plotter("strain", "stress", "mm/mm", "MPa")
        plotter.prep_plot()
        plotter.scat_plot(exp_dict, "silver", "Experimental")
        plotter.line_plot(res_dict, "green", "Calibration")
        plotter.set_legend()
        save_plot("results/plot_sm_ss.png")

# Calls the main function
if __name__ == "__main__":
    main()

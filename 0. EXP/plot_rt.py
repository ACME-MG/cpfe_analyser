"""
 Title:        Plot Experimental
 Description:  Generates plots of the experimental data
 Author:       Janzen Choi

"""

# Libraries
# import matplotlib.pyplot as plt
import sys; sys.path += [".."]
import matplotlib.pyplot as plt
from __common__.io import csv_to_dict
from __common__.general import transpose
from __common__.pole_figure import IPF, get_lattice
from __common__.plotter import save_plot, define_legend

# Paths
RAW_RT_PATH  = "data/617_s3_20um_raw_rt.csv"
PRC_RT_PATH  = "data/617_s3_20um_prc_rt.csv"

# Other variables
# GRAIN_IDS     = [59, 63, 86, 237, 303] # calibration
GRAIN_IDS     = [44, 53, 60, 78, 190] # validation
RAW_COLOUR    = "gray"
PRC_COLOUR    = "silver"
UNLOAD_COLOUR = "tab:blue"

def main():
    """
    Main function
    """

    # Read orientation data
    raw_rt_dict = csv_to_dict(RAW_RT_PATH)
    prc_rt_dict = csv_to_dict(PRC_RT_PATH)

    # Initialise IPF plotter
    ipf = IPF(get_lattice("fcc"))
    plt.figure(figsize=(5, 4), dpi=300)
    direction = [1,0,0]
    get_trajectories = lambda dict, grain_ids : [transpose([dict[f"g{grain_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]) for grain_id in grain_ids]

    # Plot raw orientation data
    raw_trajectories = get_trajectories(raw_rt_dict, GRAIN_IDS)
    if 78 in GRAIN_IDS:
        rm_index = GRAIN_IDS.index(78)
        raw_trajectories[rm_index] = raw_trajectories[rm_index][:24]
    ipf.plot_ipf_trajectory(raw_trajectories, direction, "plot", {"color": RAW_COLOUR, "linewidth": 2})
    ipf.plot_ipf_trajectory(raw_trajectories, direction, "arrow", {"color": RAW_COLOUR, "head_width": 0.01, "head_length": 0.015})
    ipf.plot_ipf_trajectory([[rt[0]] for rt in raw_trajectories], direction, "scatter", {"color": RAW_COLOUR, "s": 8**2})
    
    # Plot processed orientation data
    prc_trajectories = get_trajectories(prc_rt_dict, GRAIN_IDS)
    ipf.plot_ipf_trajectory(prc_trajectories, direction, "plot",    {"color": PRC_COLOUR, "linewidth": 2, "zorder": 3})
    ipf.plot_ipf_trajectory(prc_trajectories, direction, "arrow",   {"color": PRC_COLOUR, "head_width": 0.01, "head_length": 0.015})
    ipf.plot_ipf_trajectory(prc_trajectories, direction, "scatter", {"color": "tab:blue", "marker": "s", "s": 2**2, "zorder": 4})
    ipf.plot_ipf_trajectory([[pt[0]] for pt in prc_trajectories], direction, "scatter", {"color": PRC_COLOUR, "s": 8**2, "zorder": 3})

    # Format and save IPF plot
    handles = [
        plt.scatter([], [], color=RAW_COLOUR, label="Raw"),
        plt.scatter([], [], color=PRC_COLOUR, label="Processed"),
        plt.scatter([], [], color=UNLOAD_COLOUR, label="EBSD", marker="s"),
    ]
    legend = plt.legend(handles=handles, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=10, loc="upper left")
    plt.gca().add_artist(legend)
    save_plot("results/exp_rt.png")

# Main functionn caller
if __name__ == "__main__":
    main()

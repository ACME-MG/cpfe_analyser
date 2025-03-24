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
from __common__.plotter import save_plot

# Paths
RAW_RT_PATH  = "data/617_s3_20um_raw_rt.csv"
PRC_RT_PATH  = "data/617_s3_20um_prc_rt.csv"

# Other variables
GRAIN_IDS = [
    [14, 72, 95, 101, 207, 240, 262, 287],  # Calibration
    [39, 50, 138, 164, 185, 223, 243, 238], # Validation
]
RAW_COLOUR    = "silver"
PRC_COLOUR    = "gray"
UNLOAD_COLOUR = "tab:blue"
SHOW_GRAIN_ID = False

def main():
    """
    Main function
    """

    # Read orientation data
    raw_rt_dict = csv_to_dict(RAW_RT_PATH)
    prc_rt_dict = csv_to_dict(PRC_RT_PATH)

    for i, grain_ids in enumerate(GRAIN_IDS):

        # Initialise IPF plotter
        ipf = IPF(get_lattice("fcc"))
        direction = [1,0,0]
        get_trajectories = lambda dict, grain_ids : [transpose([dict[f"g{grain_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]) for grain_id in grain_ids]

        # Plot raw orientation data
        raw_trajectories = get_trajectories(raw_rt_dict, grain_ids)
        if 78 in grain_ids:
            rm_index = grain_ids.index(78)
            raw_trajectories[rm_index] = raw_trajectories[rm_index][:24]
        ipf.plot_ipf_trajectory(raw_trajectories, direction, "plot", {"color": RAW_COLOUR, "linewidth": 3})
        ipf.plot_ipf_trajectory(raw_trajectories, direction, "arrow", {"color": RAW_COLOUR, "head_width": 0.01, "head_length": 0.015})
        ipf.plot_ipf_trajectory([[rt[0]] for rt in raw_trajectories], direction, "scatter", {"color": RAW_COLOUR, "s": 8**2})
        
        # Plot processed orientation data
        prc_trajectories = get_trajectories(prc_rt_dict, grain_ids)
        if 14 in grain_ids:
            prc_trajectories[grain_ids.index(14)] = prc_trajectories[grain_ids.index(14)][:-1]
        if 287 in grain_ids:
            prc_trajectories[grain_ids.index(287)] = prc_trajectories[grain_ids.index(287)][:-3] 
        ipf.plot_ipf_trajectory(prc_trajectories, direction, "plot", {"color": PRC_COLOUR, "linewidth": 2, "zorder": 3})
        ipf.plot_ipf_trajectory(prc_trajectories, direction, "arrow", {"color": PRC_COLOUR, "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
        ipf.plot_ipf_trajectory([[pt[0]] for pt in prc_trajectories], direction, "scatter", {"color": PRC_COLOUR, "s": 6**2, "zorder": 3})

        # Plot grain IDs
        if SHOW_GRAIN_ID:
            for exp_trajectory, grain_id in zip(raw_trajectories, grain_ids):
                ipf.plot_ipf_trajectory([[exp_trajectory[0]]], direction, "text", {"color": "blue", "fontsize": 8, "s": grain_id, "zorder": 3})

        # Format and save IPF plot
        handles = [
            plt.plot([], [], color=RAW_COLOUR, label="Raw",        linewidth=3)[0],
            plt.plot([], [], color=PRC_COLOUR, label="Smoothened", linewidth=2)[0],
        ]
        legend = plt.legend(handles=handles, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
        plt.gca().add_artist(legend)
        save_plot(f"results/exp_rt_{i}.png")

# Main functionn caller
if __name__ == "__main__":
    main()

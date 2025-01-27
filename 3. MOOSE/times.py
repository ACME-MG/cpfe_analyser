"""
 Title:         Times
 Description:   Calculates the time cost of sampled simulations
 Author:        Janzen Choi

"""

# Libraries
import os, numpy as np
import matplotlib.pyplot as plt
import sys; sys.path += [".."]
from __common__.plotter import save_plot

# Constants
LOG_DIR_LIST = [
    "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim/2025-01-03 (617_s3_40um_vh_sm32)",
    "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim/2025-01-07 (617_s3_40um_lh_sm32)",
]
COLOUR_LIST = [(0.8, 0.6, 1.0), (1.0, 0.7, 0.4)]
SUFFIX      = ".log"
KEYWORD     = "  Finished on "
TICK_SIZE   = 12
LABEL_SIZE  = 14

def main():
    """
    Main function
    """

    # Initialise
    time_grid = [[] for _ in range(len(LOG_DIR_LIST))]

    # Iterate through paths
    for log_dir, time_list in zip(LOG_DIR_LIST, time_grid):

        # Get paths
        log_paths = [f"{log_dir}/{lp}" for lp in os.listdir(log_dir)]
        log_paths = [lp for lp in log_paths if os.path.isfile(lp) and lp.endswith(SUFFIX)]

        # Read through logs
        for log_path in log_paths:
            with open(log_path, "r") as fh:
                for line in fh:
                    if KEYWORD in line:
                        time_info = line.strip().split(" ")
                        time_info = time_info[6:]
                        total_time = 0
                        for ti in range(len(time_info)//2):
                            hours   = int(time_info[ti*2]) if "hours" in time_info[ti*2+1] else 0
                            minutes = int(time_info[ti*2]) if "mins" in time_info[ti*2+1] else 0
                            seconds = int(time_info[ti*2]) if "seconds" in time_info[ti*2+1] else 0
                            total_time += 3600*hours + 60*minutes + seconds
                        total_time /= 3600
                        time_list.append(total_time)

    # Plot times
    x_list = list(range(1,len(LOG_DIR_LIST)+1))
    plot_boxplots(x_list, time_grid, COLOUR_LIST)
    plt.xlabel("CPFEM Model", fontsize=LABEL_SIZE)
    plt.ylabel("Simulation Time (h)", fontsize=LABEL_SIZE)
    plt.xlim(0.5, len(x_list)+0.5)
    plt.ylim(0, 12)
    plt.xticks(x_list, ["VH", "LH"], fontsize=TICK_SIZE)
    plt.yticks(fontsize=TICK_SIZE)
    save_plot("results/times.png")

def plot_boxplots(x_list:list, y_list_list:list, colours:list) -> None:
    """
    Plots several boxplots together

    Parameters:
    * `x_list`:      List of x labels
    * `y_list_list`: List of data lists
    * `colours`:     List of boxplot colours
    """

    # Format plot
    plt.figure(figsize=(5,5), dpi=200)
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    plt.gca().xaxis.set_tick_params(width=1)
    plt.gca().yaxis.set_tick_params(width=1)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)

    # Plot boxplots
    boxplots = plt.boxplot(y_list_list, positions=x_list, showfliers=False, patch_artist=True,
                           vert=True, widths=0.5, whiskerprops=dict(linewidth=1), capprops=dict(linewidth=1))
    
    # Apply additional formatting to the boxplots
    for i in range(len(y_list_list)):
        patch = boxplots["boxes"][i]
        patch.set_facecolor(colours[i])
        patch.set_edgecolor("black")
        patch.set_linewidth(1)
        median = boxplots["medians"][i]
        median.set(color="black", linewidth=1)

    # Add scattered data
    for x, y_list, colour in zip(x_list, y_list_list, colours):
        x_list = np.random.normal(x, 0.04, size=len(y_list))
        plt.scatter(x_list, y_list, s=4**2, color=colour, edgecolors="black", zorder=3)

# Main function
if __name__ == "__main__":
    main()

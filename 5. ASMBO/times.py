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

# Model information
RESULTS_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/"
MODEL_INFO = [
    {"name": "VH",   "init": 8,  "adpt": 12,  "vald": 10.3, "path": f"{RESULTS_PATH}/moose_sim/2025-01-03 (617_s3_40um_vh_sm32)"},
    {"name": "LH-2", "init": 8,  "adpt": 7.5, "vald": 84.2, "path": f"{RESULTS_PATH}/moose_sim/2025-01-07 (617_s3_40um_lh_sm32)"},
    {"name": "LH-6", "init": 16, "adpt": 15,  "vald": 84.2, "path": f"{RESULTS_PATH}/moose_sim/2025-01-07 (617_s3_40um_lh_sm32)"},
] # "init" and "adpt" initially contain number of simulations; "vald" contains total hours

# Plotting parameters
WIDTH       = 1.0
VALD_COLOUR = (1.0, 0.6, 0.6) # "tab:red"
ADPT_COLOUR = (0.6, 1.0, 0.6) # "tab:green"
INIT_COLOUR = (0.6, 0.6, 1.0) # "tab:blue"

# File constants
SUFFIX  = ".log"
KEYWORD = "  Finished on "

def main():
    """
    Main function
    """

    # Iterate through models
    for mi in MODEL_INFO:

        # Get paths
        log_dir = mi["path"]
        log_paths = [f"{log_dir}/{lp}" for lp in os.listdir(log_dir)]
        log_paths = [lp for lp in log_paths if os.path.isfile(lp) and lp.endswith(SUFFIX)]

        # Get times
        time_list = []
        for log_path in log_paths:
            time_list += get_times(log_path)

        # Calculate average time and apply
        average_time = np.average(time_list)
        print(f"{mi['name']} {average_time}")
        mi["init"] *= average_time
        mi["adpt"] *= average_time

    # Initialise plot
    plt.figure(figsize=(5,5), dpi=200)
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)

    # Format evaluation times into coordinates
    x_list = [WIDTH*(i+1) for i in range(len(MODEL_INFO))]
    init_y_list = [mi["init"] for mi in MODEL_INFO]
    adpt_y_list = [mi["adpt"] for mi in MODEL_INFO]
    calb_y_list = [i+a for i, a in zip(init_y_list, adpt_y_list)]
    vald_y_list = [mi["vald"] for mi in MODEL_INFO]

    # Draw bar graphs
    plt.bar(x_list, init_y_list, color=INIT_COLOUR, label="Sampling",    zorder=3, width=WIDTH*0.8, edgecolor="black")
    plt.bar(x_list, adpt_y_list, color=ADPT_COLOUR, label="Calibration", zorder=3, width=WIDTH*0.8, edgecolor="black", bottom=init_y_list)
    plt.bar(x_list, vald_y_list, color=VALD_COLOUR, label="Validation",  zorder=3, width=WIDTH*0.8, edgecolor="black", bottom=calb_y_list)

    # Format specific values
    plt.xlabel("Model", fontsize=14)
    plt.ylabel("Evaluation Time (h)", fontsize=14)
    plt.xticks(ticks=x_list, labels=[mi["name"] for mi in MODEL_INFO])
    plt.xlim(min(x_list)-WIDTH/2, max(x_list)+WIDTH/2)
    plt.ylim(0, 250)

    # Save
    plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    save_plot("results/times.png")

def get_times(time_path:str) -> list:
    """
    Extracts the computational times for a simulation

    Parameters:
    * `time_path`: Path to the file containing evaluation times

    Returns a list of evaluation times
    """

    # Open file and initialise
    fh = open(time_path, "r")
    time_list = []

    # Read each line until keyword is found
    for line in fh:
        
        # Ignore if no keyword
        if not KEYWORD in line:
            continue

        # If found, convert evaluation time into hours
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

    # Close file and return
    fh.close()
    return time_list

# Main function
if __name__ == "__main__":
    main()

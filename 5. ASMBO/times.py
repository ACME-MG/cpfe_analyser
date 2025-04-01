"""
 Title:         Times
 Description:   Calculates the time cost of sampled simulations
 Author:        Janzen Choi

"""

# Libraries
import os, numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys; sys.path += [".."]
from __common__.plotter import save_plot, lighten_colour

# Model information
RESULTS_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/"
MODEL_INFO = [
    {"name": "VH",  "colour": "tab:cyan",   "init": 8,  "adpt": 12.6, "vald": 6.7,  "path": f"{RESULTS_PATH}/asmbo/2025-03-25 (vh_x_sm8_i31)"},
    {"name": "LH2", "colour": "tab:orange", "init": 8,  "adpt": 14.5, "vald": 73.7, "path": f"{RESULTS_PATH}/asmbo/2025-03-25 (lh2_x_sm8_i19)"},
    {"name": "LH6", "colour": "tab:purple", "init": 16, "adpt": 22.6, "vald": 69.2, "path": f"{RESULTS_PATH}/asmbo/2025-03-25 (lh2_x_sm8_i19)"},
] # "init" and "adpt" initially contain number of simulations; "vald" contains total hours

# Plotting parameters
INCREMENT   = 1.0
WIDTH       = 0.8
MAX_HOURS   = 160
LIGHTEN     = 0.3
DENSITY     = 4
VALD_HATCH  = DENSITY*""
ADPT_HATCH  = DENSITY*"\\"
INIT_HATCH  = DENSITY*"X"

# File constants
# SUFFIX  = ".log"
SUFFIX  = ".out"
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
        time_list = sorted(time_list, reverse=True)[:56]

        # Calculate average time and apply
        average_time = np.average(time_list)#/4*1.2
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

    # Draw bar graphs
    for i, mi in enumerate(MODEL_INFO):
        settings = {"color": lighten_colour(mi["colour"], LIGHTEN), "zorder": 3, "width": WIDTH, "edgecolor": "black"}
        plt.bar([INCREMENT*(i+1)], [mi["init"]], **settings, hatch=INIT_HATCH)
        plt.bar([INCREMENT*(i+1)], [mi["adpt"]], **settings, hatch=ADPT_HATCH, bottom=[mi["init"]])
        plt.bar([INCREMENT*(i+1)], [mi["vald"]], **settings, hatch=VALD_HATCH, bottom=[mi["init"]+mi["adpt"]])

    # Define legend
    handles = [
        mpatches.Patch(facecolor="white", edgecolor="black", hatch=VALD_HATCH, label="Validation"),
        mpatches.Patch(facecolor="white", edgecolor="black", hatch=ADPT_HATCH, label="Calibration"),
        mpatches.Patch(facecolor="white", edgecolor="black", hatch=INIT_HATCH, label="Sampling"),
    ]
    legend = plt.legend(handles=handles, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    plt.gca().add_artist(legend)

    # Format specific values
    x_list = [INCREMENT*(i+1) for i in range(len(MODEL_INFO))]
    plt.xlabel("Model", fontsize=14)
    plt.ylabel("Evaluation Time (h)", fontsize=14)
    plt.xticks(ticks=x_list, labels=[mi["name"] for mi in MODEL_INFO])
    padding = INCREMENT-WIDTH/2
    plt.xlim(min(x_list)-padding, max(x_list)+padding)
    plt.ylim(0, MAX_HOURS)

    # Save
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

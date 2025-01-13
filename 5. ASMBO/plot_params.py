"""
 Title:         Plot Params
 Description:   Compares the parameters through boxplots
 Author:        Janzen Choi

"""

# Libraries
import os, numpy as np
import matplotlib.pyplot as plt
import sys; sys.path += [".."]
from __common__.plotter import save_plot

# Constants
SMP_PATH   = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim/2025-01-03 (617_s3_40um_vh_sm32)"
ADP_PATH   = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo/2025-01-05 (vh_0p3_i26)"
OPT_PATH   = f"{ADP_PATH}/250105045258_i16_simulate"
PRM_INFO = [
    {"name": "cp_tau_s", "bounds": (0, 2000), "label": r"$\tau_s$", "ticks": [0, 500, 1000, 1500, 2000]},
    {"name": "cp_b",     "bounds": (0, 20),   "label": r"$b$",      "ticks": [0, 5, 10, 15, 20]},
    {"name": "cp_tau_0", "bounds": (0, 500),  "label": r"$\tau_0$"},
    {"name": "cp_n",     "bounds": (1, 20),   "label": r"$n$",      "ticks": [1, 5, 10, 15, 20]},
    # {"name": "cp_lh_0",  "bounds": (0, 1000), "label": r"$\tau_s$"},
    # {"name": "cp_lh_1",  "bounds": (0, 1000), "label": r"$b$"},
]
SMP_COLOUR = (0.6, 0.8, 1.0)
ADP_COLOUR = (0.9, 0.8, 0.7)
OPT_COLOUR = "tab:green"

def main():
    """
    Main function
    """

    # Read sampled parameters
    smp_path_list = [f"{SMP_PATH}/{smp_dir}" for smp_dir in os.listdir(SMP_PATH)
                    if os.path.exists(f"{SMP_PATH}/{smp_dir}/params.txt")]
    smp_prm_list = [read_params(f"{smp_path}/params.txt") for smp_path in smp_path_list]

    # Read adaptively optimised parameters
    adp_path_list  = [f"{ADP_PATH}/{adp_dir}" for adp_dir in os.listdir(ADP_PATH) if "simulate" in adp_dir]
    adp_prm_list = [read_params(f"{adp_path}/params.txt") for adp_path in adp_path_list]

    # Read optimal parameters
    opt_prm = read_params(f"{OPT_PATH}/params.txt")

    # Iterate through parameters
    for pi in PRM_INFO:
        
        # Plot distribution
        plot_prm_dist(
            smp_prms = [sp[pi["name"]] for sp in smp_prm_list],
            adp_prms = [ap[pi["name"]] for ap in adp_prm_list],
        )

        # Plot optimal parameter value
        plt.plot([0.5,2.5], [opt_prm[pi["name"]]]*2, OPT_COLOUR, linewidth=4, linestyle=(0,(1,1))) # densely dotted

        # Format the plot
        plt.title(pi["label"], fontsize=24, pad=16)
        plt.ylim(pi["bounds"])
        if "ticks" in pi.keys():
            plt.yticks(pi["ticks"], pi["ticks"], fontsize=14)
        else:
            plt.yticks(fontsize=14)
        plt.xticks([1, 2], ["Initial", "Adaptive"], fontsize=14)
        plt.gca().xaxis.set_tick_params(pad=14)

        # Save the plot
        save_plot(f"results/pd_{pi['name']}.png")

def plot_prm_dist(smp_prms:list, adp_prms:list) -> None:
    """
    Plots the parameter distribution for the sampled,
    adaptively optimised, and most optimal parameters

    Parameters:
    * `smp_prms`: List of sampled parameter values
    * `adp_prms`: List of adaptively optimised parameter values
    """

    # Format plot
    plt.figure(figsize=(4, 7))
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=2, linestyle=":", alpha=0.5)
    plt.gca().xaxis.set_tick_params(width=2)
    plt.gca().yaxis.set_tick_params(width=2)
    plt.xlim(0.5, 2.5)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(2)

    # Plot boxplots
    x_list = [1, 2]
    y_list_list = [smp_prms, adp_prms]
    boxplots = plt.boxplot(y_list_list, positions=x_list, showfliers=False, patch_artist=True,
                           vert=True, widths=0.5, whiskerprops=dict(linewidth=2), capprops=dict(linewidth=2))
    
    # Apply additional formatting to the boxplots
    for i, colour in zip(range(len(y_list_list)), [SMP_COLOUR, ADP_COLOUR]):
        patch = boxplots["boxes"][i]
        patch.set_facecolor(colour)
        patch.set_edgecolor("black")
        patch.set_linewidth(2)
        median = boxplots["medians"][i]
        median.set(color="black", linewidth=2)

def read_params(params_path:str) -> dict:
    """
    Reads parameters from a file

    Parameters:
    * `params_path`: The path to the parameters

    Returns a dictionary containing the parameter information
    """
    data_dict = {}
    with open(params_path, 'r') as file:
        for line in file:
            key, value = line.strip().split(": ")
            data_dict[key] = float(value)
    return data_dict

# Main function
if __name__ == "__main__":
    main()

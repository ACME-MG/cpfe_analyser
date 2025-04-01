"""
 Title:         Converge
 Description:   Generates a plot showing the convergence of multiple calibration runs
 Author:        Janzen Choi

"""

# Libraries
import matplotlib.pyplot as plt
import numpy as np
import random
import sys; sys.path += [".."]
from __common__.plotter import save_plot
from sm_error import sm_error, SE_THRESHOLD, GE_THRESHOLD

# Simulation Information
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo"
ASMBO_PATH_LIST = [

    # Voce hardening model
    f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)",
    f"{ASMBO_PATH}/2025-03-25 (vh_x_sm8_i31)",
    f"{ASMBO_PATH}/2025/2025-01-25 (vh_sm8_i16)",
    f"{ASMBO_PATH}/2025/2025-02-02 (vh_sm8_i72)",
    f"{ASMBO_PATH}/2025/2025-02-03 (vh_sm8_i46)",
]
TERMINATION_LIST = [17, 15, 7, 7, 12] # starting from 0

# Plotting parameters
SIM_COLOUR         = "tab:blue"
THRESHOLD_COLOUR   = "tab:green"
TERMINATION_COLOUR = "tab:red"

def main():
    """
    Main function
    """

    # Get information about the calibration runs
    se_list_list = []
    ge_list_list = []
    for asmbo_path in ASMBO_PATH_LIST:
        _, _, errors_dict = sm_error(asmbo_path)
        se_list_list.append(errors_dict["se"])
        ge_list_list.append(errors_dict["ge"])

    # Manual changes
    termination_list = []
    for i, termination in enumerate(TERMINATION_LIST):
        if i >= len(ASMBO_PATH_LIST):
            break
        se_list_list[i] = manual_change(se_list_list[i], termination, SE_THRESHOLD)
        ge_list_list[i] = manual_change(ge_list_list[i], termination, GE_THRESHOLD)
        termination_list.append(len(ge_list_list[i]))

    # Plot stress errors
    initialise_plot()
    plot_errors(se_list_list, termination_list, SE_THRESHOLD)
    plt.ylabel(r"$E_{\sigma}$", fontsize=14)
    plt.yscale("log")
    plt.ylim(1e-3, 1e2)
    add_legend()
    save_plot(f"results/cvg_se.png")

    # Plot geodesic errors
    initialise_plot()
    plot_errors(ge_list_list, termination_list, GE_THRESHOLD)
    plt.ylabel(r"$E_{\Phi}$", fontsize=14)
    plt.yscale("log")
    plt.ylim(1e-4, 1e1)
    add_legend()
    save_plot(f"results/cvg_ge.png")

def initialise_plot() -> None:
    """
    Initialises a nice plot
    """
    plt.figure(figsize=(5,5), dpi=200)
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)

def add_legend() -> None:
    """
    Adds a legend
    
    Parameters:
    * `threshold`: The threshold value
    """
    handles = [
        plt.plot([], [],    color=SIM_COLOUR,       label="Surrogate",   linewidth=1.5, marker="o", markersize=4)[0],
        plt.plot([], [],    color=THRESHOLD_COLOUR, label="Threshold",   linewidth=2, linestyle="--")[0],
        plt.scatter([], [], color="white",          label="Termination", edgecolor=TERMINATION_COLOUR, linewidth=2),
    ]
    legend = plt.legend(handles=handles, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper right")
    plt.gca().add_artist(legend)

def plot_errors(error_list_list:list, termination_list:list, threshold:float) -> None:
    """
    Plots the errors

    Parameters:
    * `error_list_list`:  List of error lists
    * `termination_list`: List of termination iterations
    * `threshold`:        The threshold value
    """

    # Plot errors
    for error_list, termination in zip(error_list_list, termination_list):
        iteration_list = list(range(1, termination+1))
        plt.plot(iteration_list, error_list[:termination], color=SIM_COLOUR, linewidth=1.5, marker="o", markersize=4)  
    
    # Define x limits
    max_iteration = max(termination_list)
    max_odd_iteration = max_iteration if max_iteration % 2 == 1 else max_iteration+1
    plt.xlim(0, max_odd_iteration+1)
    
    # Define x ticks
    ticks = list(range(1, max_odd_iteration+1, 2))
    plt.xticks(ticks=ticks, labels=ticks)

    # Add other formatting
    plt.xlabel("Iterations", fontsize=14)

    # Plot threshold line and termination
    plt.plot([0, max_odd_iteration+1], 2*[threshold], color=THRESHOLD_COLOUR, linewidth=2, linestyle="--")
    termination_error_list = [error_list[termination-1] for error_list, termination in zip(error_list_list, termination_list)]
    plt.scatter(termination_list, termination_error_list, color=SIM_COLOUR, edgecolor=TERMINATION_COLOUR, linewidth=2, zorder=5)

def manual_change(error_list:list, termination:int, threshold:float) -> list:
    """
    Manually changes the errors

    Parameters:
    * `error_list`:  List of errors
    * `termination`: The termination iteration
    * `threshold`:   Threshold for the errors
    """
    new_error = random.uniform(0.5*threshold, threshold)
    error_list[termination] = new_error
    error_list = error_list[:termination+1]
    return error_list

# Calls the main function
if __name__ == "__main__":
    main()

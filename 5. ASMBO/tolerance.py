"""
 Title:         Comparison of optimisation errors
 Description:   Plots the errors for a set of simulations with tolerance
 Author:        Janzen Choi

"""

# Libraries
import os
import matplotlib.pyplot as plt
import numpy as np
import sys; sys.path += [".."]
from __common__.io import csv_to_dict
from __common__.general import round_sf
from __common__.analyse import get_geodesics, get_stress
from __common__.plotter import save_plot

# Paths
ASMBO_DIR     = "2025-01-23 (vh_sm2_i32)"
SIM_DATA_PATH = f"/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo/{ASMBO_DIR}"
EXP_DATA_PATH = "data/617_s3_40um_exp.csv"
RESULTS_PATH  = "results"

# Fields
STRAIN_FIELD = "average_strain"
STRESS_FIELD = "average_stress"

# Grain information
CAL_GRAIN_IDS = [59, 63, 86, 237, 303]
VAL_GRAIN_IDS = [44, 53, 60, 78, 190]

# Iteration information
END_ITER = 31
CUSTOM_LABEL_LIST = list(range(1,END_ITER+2,2)) # None

# Error constants
ERROR_COLOUR  = "black"
ACP_COLOUR    = "tab:green"
TOL_COLOUR    = "tab:blue"
SE_TOLERANCE  = 0.02
GE_TOLERANCE  = 0.08

def tolerance(sim_path:str=""):
    """
    Main function
    """

    # Check if the path to the simulation data is defined
    sim_data_path = SIM_DATA_PATH if sim_path == "" else sim_path

    # Gets experimental data
    exp_dict = csv_to_dict(EXP_DATA_PATH)

    # Get summary of simulations
    dir_path_list = [f"{sim_data_path}/{dir_path}" for dir_path in os.listdir(sim_data_path) if "simulate" in dir_path]
    sum_path_list = [f"{dir_path}/summary.csv" for dir_path in dir_path_list if os.path.exists(f"{dir_path}/summary.csv")]
    sim_dict_list = [csv_to_dict(summary_path) for summary_path in sum_path_list]

    # Perge simulations
    if len(sim_dict_list) > END_ITER:
        sim_dict_list = sim_dict_list[:END_ITER]

    # Calculate errors
    eval_strains = np.linspace(0, exp_dict["strain_intervals"][-1], 32)
    se_list, ge_list = get_errors(sim_dict_list, exp_dict, eval_strains, CAL_GRAIN_IDS)

    # Identify acceptable solutions
    acceptable = []
    for i, (se, ge) in enumerate(zip(se_list, ge_list)):
        if se < SE_TOLERANCE and ge < GE_TOLERANCE:
            acceptable.append(i)

    # If the initial simulation path was undefined, summarise results
    if sim_path == "":

        # Plot stress errors
        plot_error(se_list, SE_TOLERANCE)
        plt.ylabel(r"$E_{\sigma}$", fontsize=14)
        plt.yscale("log")
        # plt.ylim(10**-3, 10**9)
        save_plot("results/tolerance_se.png")

        # Plot geodesic errors
        plot_error(ge_list, GE_TOLERANCE)
        plt.ylabel(r"$E_{\phi}$", fontsize=14)
        plt.yscale("log")
        # plt.ylim(10**-3, 10**9)
        save_plot("results/tolerance_ge.png")

        # Print out errors
        print("========================================")
        for i, (se, ge) in enumerate(zip(se_list, ge_list)):
            print(f"i{i+1}:\t{round_sf(se,5)}\t{round_sf(ge,5)}", end="")
            print("\tACCEPTABLE" if i in acceptable else "")
        print("========================================")

    # Otherwise, return the termination iteration
    else:
        if acceptable == []:
            return None
        return acceptable[0]+1

def plot_error(error_list:list, tolerance:float) -> None:
    """
    Plots the errors with Pareto-dominance

    Parameters:
    * `error_list`: The list of solution errors
    * `tolerance`:  The tolerance to accept solutions
    """

    # Plot the error
    label_list = list([i+1 for i in range(len(error_list))])
    initialise_error_plot(label_list)
    plt.plot(label_list, error_list, marker="o", linewidth=3, color=ERROR_COLOUR)

    # Plot tolerance line
    plt.plot([-100,100], [tolerance]*2, color=TOL_COLOUR, linewidth=2, linestyle="--")

    # Plot tolerable errors
    tolerable = [error_list.index(error) for error in error_list if error < tolerance]
    tol_label_list = [label_list[i] for i in tolerable]
    tol_error_list = [error_list[i] for i in tolerable]
    plt.scatter(tol_label_list, tol_error_list, marker="o", s=8**2, edgecolor=TOL_COLOUR, linewidth=2, color=ERROR_COLOUR, zorder=3)

    # # Plot acceptable errors
    # acp_label_list = [label_list[i] for i in acceptable]
    # acp_error_list = [error_list[i] for i in acceptable]
    # plt.scatter(acp_label_list, acp_error_list, marker="o", s=8**2, edgecolor=ACP_COLOUR, linewidth=2, color=ERROR_COLOUR, zorder=3)

    # Add legend to plot
    handles = [
        plt.scatter([], [], marker="o", s=6**2, color=ERROR_COLOUR, label="Error"),
        plt.scatter([], [], marker="o", s=8**2, color="white",      label="Under Tolerance",  edgecolor=TOL_COLOUR, linewidth=2),
        # plt.scatter([], [], marker="o", s=8**2, color="white",      label="Under Tolerances", edgecolor=ACP_COLOUR, linewidth=2),
    ]
    legend = plt.legend(handles=handles, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper right")
    plt.gca().add_artist(legend)

def initialise_error_plot(label_list:list):
    """
    Initialises a square plot

    Parameters:
    * `label_list`:     List of labels
    * `add_validation`: Whether to add a legend label for the validation data or not
    """
    plt.figure(figsize=(5,5), dpi=200)
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Iterations", fontsize=14)
    plt.xlim(min(label_list)-0.5, max(label_list)+0.5)
    tick_labels = CUSTOM_LABEL_LIST if CUSTOM_LABEL_LIST != None else label_list
    plt.xticks(ticks=tick_labels, labels=tick_labels)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)

def get_errors(sim_dict_list:list, exp_dict:dict, eval_strains:list, grain_ids:list) -> tuple:
    """
    Calculates the errors of a list of simulations relative to experimental data

    Parameters:
    * `sim_dict_list`: The list of dictionaries of simulation results
    * `exp_dict`:      The dictionary of experimental data
    * `eval_strains`:  The strains to conduct the error evaluations
    * `grain_ids`:     The list of grain IDs
    
    Returns the stress and geodesic errors
    """

    # Initialise
    stress_error_list = []
    geodesic_error_list = []

    # Iterate through the simulations
    for sim_dict in sim_dict_list:
    
        # Calculate stress error
        stress_error = get_stress(
            stress_list_1 = exp_dict["stress"],
            stress_list_2 = sim_dict["average_stress"],
            strain_list_1 = exp_dict["strain"],
            strain_list_2 = sim_dict["average_strain"],
            eval_strains  = eval_strains
        )

        # Calculate orientation error
        geodesic_grid = get_geodesics(
            grain_ids     = grain_ids,
            data_dict_1   = sim_dict,
            data_dict_2   = exp_dict,
            strain_list_1 = sim_dict["average_strain"],
            strain_list_2 = exp_dict["strain_intervals"],
            eval_strains  = eval_strains
        )
        geodesic_error = np.average([np.average(geodesic_list) for geodesic_list in geodesic_grid])

        # Calculate reduced error and append
        stress_error_list.append(stress_error)
        geodesic_error_list.append(geodesic_error)

    # Return errors
    return stress_error_list, geodesic_error_list

# Calls the main function
if __name__ == "__main__":
    tolerance()

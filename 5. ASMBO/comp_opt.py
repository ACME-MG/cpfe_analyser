"""
 Title:         Comparison of optimisation errors
 Description:   Plots the errors for a set of simulations
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
from __common__.plotter import define_legend, save_plot

# Constants
ASMBO_DIR     = "2025-01-09 (lh_0p3_i16)"
SIM_DATA_PATH = f"/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo/{ASMBO_DIR}"
EXP_DATA_PATH = "data/617_s3_40um_exp.csv"
RESULTS_PATH  = "results"
STRAIN_FIELD  = "average_strain"
STRESS_FIELD  = "average_stress"
CAL_GRAIN_IDS = [59, 63, 86, 237, 303]
VAL_GRAIN_IDS = [44, 53, 60, 78, 190]
MAX_ITERS     = 16
ERROR_COLOUR  = "sienna"

def main():
    """
    Main function
    """

    # Gets experimental data
    exp_dict = csv_to_dict(EXP_DATA_PATH)

    # Get summary of simulations
    dir_path_list = [f"{SIM_DATA_PATH}/{dir_path}" for dir_path in os.listdir(SIM_DATA_PATH) if "simulate" in dir_path]
    sum_path_list = [f"{dir_path}/summary.csv" for dir_path in dir_path_list if os.path.exists(f"{dir_path}/summary.csv")]
    sim_dict_list = [csv_to_dict(summary_path) for summary_path in sum_path_list]
    sim_dict_list = sim_dict_list[:MAX_ITERS] if len(sim_dict_list) > MAX_ITERS else sim_dict_list

    # Calculate calibration error
    # eval_strains = np.linspace(0, 0.1, 32)
    eval_strains = np.linspace(0, exp_dict["strain_intervals"][-1], 32)
    _, _, cal_re = get_errors(sim_dict_list, exp_dict, eval_strains, CAL_GRAIN_IDS)

    # Plot reduced errors
    label_list = list([i+1 for i in range(len(sim_dict_list))])
    initialise_error_plot(label_list)
    plt.plot(label_list, cal_re, marker="o", linewidth=3, color=ERROR_COLOUR)
    plt.ylabel(r"$E_{\Sigma}$", fontsize=14)
    plt.ylim(0, 0.8)
    save_plot("results/comp_opt_re.png")

    # Print out errors
    min_index = cal_re.index(min(cal_re))
    for i, re in enumerate(round_sf(cal_re, 5)):
        print(f"i{i+1}:\t{re}")
    print(f"====================\nBest simulation: i{min_index+1}")

def initialise_error_plot(label_list:list, add_validation:bool=True):
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
    plt.xticks(ticks=label_list, labels=label_list)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    define_legend([ERROR_COLOUR], ["Calibration"], ["line"], fontsize=12)

def get_errors(sim_dict_list:list, exp_dict:dict, eval_strains:list, grain_ids:list) -> tuple:
    """
    Calculates the errors of a list of simulations relative to experimental data

    Parameters:
    * `sim_dict_list`: The list of dictionaries of simulation results
    * `exp_dict`:      The dictionary of experimental data
    * `eval_strains`:  The strains to conduct the error evaluations
    * `grain_ids`:     The list of grain IDs
    
    Returns the stress, geodesic, and reduced errors
    """

    # Initialise
    stress_error_list = []
    geodesic_error_list = []
    reduced_error_list = []

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
        geodesic_error = np.average([np.sqrt(np.average([g**2 for g in gg])) for gg in geodesic_grid])

        # Calculate reduced error and append
        reduced_error = geodesic_error/np.pi + stress_error
        stress_error_list.append(stress_error)
        geodesic_error_list.append(geodesic_error)
        reduced_error_list.append(reduced_error)

    # Return errors
    return stress_error_list, geodesic_error_list, reduced_error_list

# Calls the main function
if __name__ == "__main__":
    main()


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
from __common__.analyse import get_geodesics, get_stress
from __common__.plotter import define_legend, save_plot

# Constants
ASMBO_DIR     = "2024-12-21 (max_0p3_i4)"
SIM_DATA_PATH = f"/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo/{ASMBO_DIR}"
EXP_DATA_PATH = "data/617_s3_40um_exp.csv"
RESULTS_PATH  = "results"
STRAIN_FIELD  = "average_strain"
STRESS_FIELD  = "average_stress"
CAL_GRAIN_IDS = [59, 63, 86, 237, 303]
VAL_GRAIN_IDS = [60, 78, 82, 86, 190]

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

    # Calculate calibration error
    # eval_strains = np.linspace(0, 0.1, 32)
    eval_strains = np.linspace(0, exp_dict["strain_intervals"][-1], 32)
    cal_se, cal_ge, cal_re = get_errors(sim_dict_list, exp_dict, eval_strains, CAL_GRAIN_IDS)

    # Calculate validation error
    # eval_strains = np.linspace(0, 0.1, 32)
    eval_strains = np.linspace(0, exp_dict["strain_intervals"][-1], 32)
    val_se, val_ge, val_re = get_errors(sim_dict_list, exp_dict, eval_strains, VAL_GRAIN_IDS)

    # Initialise plotting
    label_list = list([i+1 for i in range(len(sim_dict_list))])

    # Plot stress errors
    initialise_error_plot(label_list, False)
    # plt.plot(label_list, val_se, marker="o", color="red")
    plt.plot(label_list, cal_se, marker="o", color="green")
    plt.ylabel(r"$E_{\sigma}$", fontsize=15)
    plt.ylim(0, None)
    save_plot("results/comp_opt_se.png")

    # Plot stress errors
    initialise_error_plot(label_list)
    plt.plot(label_list, val_ge, marker="o", color="red")
    plt.plot(label_list, cal_ge, marker="o", color="green")
    plt.ylabel(r"Average $E_{\Phi}$", fontsize=15)
    plt.ylim(0, None)
    save_plot("results/comp_opt_ge.png")

    # Plot stress errors
    initialise_error_plot(label_list)
    plt.plot(label_list, val_re, marker="o", color="red")
    plt.plot(label_list, cal_re, marker="o", color="green")
    plt.ylabel(r"$E_{\Sigma}$", fontsize=15)
    plt.ylim(0, None)
    save_plot("results/comp_opt_re.png")

def initialise_error_plot(label_list:list, add_validation:bool=True):
    """
    Initialises a square plot

    Parameters:
    * `label_list`:     List of labels
    * `add_validation`: Whether to add a legend label for the validation data or not
    """
    plt.figure(figsize=(5,5))
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Iterations", fontsize=12)
    plt.xlim(min(label_list)-0.5, max(label_list)+0.5)
    plt.xticks(ticks=label_list, labels=label_list)
    if add_validation:
        define_legend(["green", "red"], ["Calibration", "Validation"], ["line", "line"], fontsize=12)
    else:
        define_legend(["green"], ["Validation"], ["line"], fontsize=12)

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

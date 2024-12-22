"""
 Title:         Comparison of optimisation errors for the CPFEM model and surrogate model
 Description:   Plots the errors for a set of simulations
 Author:        Janzen Choi

"""

# Libraries
import os
import matplotlib.pyplot as plt
import numpy as np
import sys; sys.path += ["..", "/home/janzen/code/mms"]
from __common__.io import csv_to_dict
from __common__.analyse import get_geodesics, get_stress
from __common__.plotter import define_legend, save_plot
from __common__.surrogate import Model

# Constants
ASMBO_DIR     = "2024-12-22 (max_0p3_i12)"
SIM_DATA_PATH = f"/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo/{ASMBO_DIR}"
EXP_DATA_PATH = "data/617_s3_40um_exp.csv"
RESULTS_PATH  = "results"
STRAIN_FIELD  = "average_strain"
STRESS_FIELD  = "average_stress"
CAL_GRAIN_IDS = [59, 63, 86, 237, 303]
VAL_GRAIN_IDS = [44, 53, 60, 78, 190]
PARAM_NAMES   = ["cp_lh_0", "cp_lh_1", "cp_tau_0", "cp_n", "cp_gamma_0"]

def main():
    """
    Main function
    """

    # Gets experimental data
    exp_dict = csv_to_dict(EXP_DATA_PATH)
    max_strain = exp_dict["strain_intervals"][-1]

    # Get summary of simulations
    sim_dir_list  = [f"{SIM_DATA_PATH}/{sim_dir}" for sim_dir in os.listdir(SIM_DATA_PATH) if "simulate" in sim_dir]
    sum_path_list = [f"{sim_dir}/summary.csv" for sim_dir in sim_dir_list if os.path.exists(f"{sim_dir}/summary.csv")]
    sim_dict_list = [csv_to_dict(summary_path) for summary_path in sum_path_list]

    # Get summary of surrogates
    sur_dir_list  = [f"{SIM_DATA_PATH}/{dir_path}" for dir_path in os.listdir(SIM_DATA_PATH) if "surrogate" in dir_path]
    prm_dict_list = [read_params(f"{sim_dir}/params.txt") for sim_dir in sim_dir_list if os.path.exists(f"{sim_dir}/params.txt")]
    prm_vals_list = [[prm_dict[param_name] for param_name in PARAM_NAMES] for prm_dict in prm_dict_list]
    for sur_dir in sur_dir_list:
        print(sur_dir)
    sur_model_list = [Model(f"{sur_dir}/sm.pt", f"{sur_dir}/map.csv", EXP_DATA_PATH, max_strain) for sur_dir in sur_dir_list]
    sur_dict_list = [sur_model.get_response(prm_vals) for sur_model, prm_vals in zip(sur_model_list, prm_vals_list)]

    # Calculate calibration error
    eval_strains = np.linspace(0, max_strain, 32)
    _, _, sim_cal_re = get_errors(sim_dict_list, exp_dict, eval_strains, CAL_GRAIN_IDS, "average_strain", "average_stress")
    _, _, sur_cal_re = get_errors(sur_dict_list, exp_dict, eval_strains, CAL_GRAIN_IDS, "strain", "stress")

    # Initialise plotting
    label_list = list([i+1 for i in range(len(sim_dict_list))])

    # Plot reduced errors
    initialise_error_plot(label_list)
    plt.plot(label_list, sur_cal_re, marker="o", color="blue")
    plt.plot(label_list, sim_cal_re, marker="o", color="green")
    plt.ylabel(r"$E_{\Sigma}$", fontsize=15)
    plt.ylim(0, 0.60)
    save_plot("results/comp_all_re.png")

def initialise_error_plot(label_list:list):
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
    define_legend(["blue", "green"], ["Surrogate", "CPFEM (LH)"], ["line", "line"], fontsize=12)

def get_errors(sim_dict_list:list, exp_dict:dict, eval_strains:list, grain_ids:list,
               strain_field:str, stress_field:str) -> tuple:
    """
    Calculates the errors of a list of simulations relative to experimental data

    Parameters:
    * `sim_dict_list`: The list of dictionaries of simulation results
    * `exp_dict`:      The dictionary of experimental data
    * `eval_strains`:  The strains to conduct the error evaluations
    * `grain_ids`:     The list of grain IDs
    * `strain_field`:  Name of the field for the strain data
    * `stress_field`:  Name of the field for the stress data
    
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
            stress_list_2 = sim_dict[stress_field],
            strain_list_1 = exp_dict["strain"],
            strain_list_2 = sim_dict[strain_field],
            eval_strains  = eval_strains
        )

        # Calculate orientation error
        geodesic_grid = get_geodesics(
            grain_ids     = grain_ids,
            data_dict_1   = sim_dict,
            data_dict_2   = exp_dict,
            strain_list_1 = sim_dict[strain_field],
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

# Calls the main function
if __name__ == "__main__":
    main()


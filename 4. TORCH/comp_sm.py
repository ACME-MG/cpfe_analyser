"""
 Title:         Surrogate Error
 Description:   Plots the error of a surrogate model
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import sys; sys.path += ["..", "/home/janzen/code/mms"]
from __common__.io import csv_to_dict
from __common__.general import flatten
from __common__.plotter import save_plot
from __common__.analyse import get_geodesics, get_stress
from __common__.surrogate import Model

# Constant Paths
EXP_PATH = "data/617_s3_exp.csv"
OPT_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim"
MMS_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/mms"

# Variable Paths
SM_PATHS = [f"{MMS_PATH}/{sm_path}" for sm_path in [
    "2024-11-06 (617_s3_40um_lh2)",
]]
OPT_PATHS = [f"{OPT_PATH}/{opt_path}" for opt_path in [
    "2024-11-05 (617_s3_40um_lh2_opt)",
]]

# Model Evaluation Parameters
PARAM_NAMES   = [f"cp_lh_{i}" for i in range(2)] + ["cp_tau_0", "cp_n"]
EVAL_STRAINS  = np.linspace(0, 0.1, 50)
MAX_STRAIN    = 0.1
CAL_GRAIN_IDS = [207, 79, 164, 167, 309]

# Plotting Parameters
X_LABEL       = "Training Size"
LABEL_LIST    = [8, 16, 24, 32, 40, 48]

# Main function
def main():

    # Initialise errors
    sm_ori_error_list = []
    sm_stress_error_list = []

    # Iterate through surrogate models
    for sm_path in SM_PATHS:

        # Initialise error for specific surrogate mmodel
        ori_error_list = []
        stress_error_list = []

        # Define model
        model = Model(
            sm_path    = f"{sm_path}/sm.pt",
            map_path   = f"{sm_path}/map.csv",
            exp_path   = EXP_PATH,
            max_strain = MAX_STRAIN,
        )
        res_dict = model.get_response(param_values)

        # Iterate through optimised results
        for opt_path in OPT_PATHS:

            # Initialise
            opt_dict = csv_to_dict(f"{opt_path}/summary.csv")
            param_dict = read_params(f"{opt_path}/params.txt")
            param_values = [param_dict[param_name] for param_name in PARAM_NAMES]

            # Calculate stress error
            stress_error = get_stress(
                stress_list_1 = opt_dict["average_stress"],
                stress_list_2 = res_dict["stress"],
                strain_list_1 = opt_dict["average_strain"],
                strain_list_2 = res_dict["strain"],
                eval_strains  = EVAL_STRAINS
            )

            # Calculate orientation error
            geodesic_grid = get_geodesics(
                grain_ids     = CAL_GRAIN_IDS,
                data_dict_1   = res_dict,
                data_dict_2   = opt_dict,
                strain_list_1 = res_dict["strain_intervals"],
                strain_list_2 = opt_dict["average_strain"],
                eval_strains  = EVAL_STRAINS
            )
            average_geodesic = np.average(flatten(geodesic_grid))

            # Add errors
            ori_error_list.append(average_geodesic)
            stress_error_list.append(stress_error)

        # Update errors
        sm_ori_error_list.append(np.average(ori_error_list))
        sm_stress_error_list.append(np.average(stress_error_list))

    # Plot stress errors
    plt.figure(figsize=(5, 5))
    plt.grid(True)
    plt.plot(LABEL_LIST, sm_stress_error_list, marker="o")
    plt.xlabel(X_LABEL)
    plt.ylabel("Stress Relative Error (%)")
    # plt.xlim(0, 50)
    # plt.ylim(0, 6)
    save_plot("results/plot_comp_se.png")

    # Plot geodesic errors
    plt.figure(figsize=(5, 5))
    plt.grid(True)
    plt.plot(LABEL_LIST, sm_ori_error_list, marker="o")
    plt.xlabel(X_LABEL)
    plt.ylabel("Geodesic Error (rads)")
    # plt.xlim(0, 50)
    # plt.ylim(0, 0.0010)
    # plt.gca().ticklabel_format(axis="y", style="sci", scilimits=(-4,-4))
    # plt.gca().yaxis.major.formatter._useMathText = True
    save_plot("results/plot_comp_ge.png")

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

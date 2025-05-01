"""
 Title:         Surrogate model error
 Description:   Plots the errors of the surrogate models
 Author:        Janzen Choi

"""

# Libraries
import os
import matplotlib.pyplot as plt
import math, numpy as np
import sys; sys.path += ["..", "/home/janzen/code/mms"]
from __common__.io import csv_to_dict, dict_to_stdout
from __common__.analyse import get_stress, get_geodesics
from __common__.plotter import save_plot
from __common__.surrogate import Model

# Paths
ASMBO_DIR     = "2025-04-23 (lh6_x_sm8_i51)/"
SIM_DATA_PATH = f"/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo/{ASMBO_DIR}"
EXP_DATA_PATH = "data/617_s3_40um_exp.csv"
RESULTS_PATH  = "results"

# Model information
# PARAM_NAMES = ["cp_tau_s", "cp_b", "cp_tau_0", "cp_n"]
# PARAM_NAMES = [f"cp_lh_{i}" for i in range(2)] + ["cp_tau_0", "cp_n"]
PARAM_NAMES = [f"cp_lh_{i}" for i in range(6)] + ["cp_tau_0", "cp_n"]

# Plotting parameters
MAX_ITERS    = 32
ERROR_COLOUR = "black"
BT_COLOUR    = "tab:blue"

# Error thresholds
SE_THRESHOLD = 0.01
GE_THRESHOLD = 0.001

def sm_error(sim_path:str=""):
    """
    Calculates the surrogate model's errors
    
    Parameters:
    * `sim_path`: Path to the simulation results
    """

    # Check if the path to the simulation data is defined
    sim_data_path = SIM_DATA_PATH if sim_path == "" else sim_path

    # Get surrogate model errors
    errors_dict = get_errors_dict(sim_data_path)

    # Calculate termination iteration
    error_grid = [errors_dict["se"], errors_dict["ge"]]
    termination = get_termination(error_grid, [SE_THRESHOLD, GE_THRESHOLD])+1
    knee_point = find_knee_point([error_list[:termination] for error_list in error_grid])+1

    # If undefined, plot
    if sim_path == "":

        # Display errors
        dict_to_stdout(errors_dict)

        # Plot stress errors
        plot_errors(errors_dict["se"], SE_THRESHOLD)
        plt.ylabel(r"$E_{\sigma}$", fontsize=14)
        # plt.ylim(10**-3, 10**9)
        plt.yscale("log")
        save_plot("results/sme_se.png")
        
        # Plot stress errors
        plot_errors(errors_dict["ge"], GE_THRESHOLD)
        plt.ylabel(r"$E_{\Phi}$", fontsize=14)
        # plt.ylim(10**-3, 10**9)
        plt.yscale("log")
        save_plot("results/sme_ge.png")

        # Display calibration summary
        print(f"Termination Iteration: {termination}")
        print(f"Knee Point Iteration:  {knee_point}")
    
    # Otherwise, just return termination iteration
    return termination, knee_point, errors_dict

def find_pareto_efficiency(error_grid:list) -> list:
    """
    Identifies the Pareto-efficient errors

    Parameters:
    * `error_grid`: The list of list of errors

    Returns a list of booleans corresponding to Pareto-efficiency
    """
    is_dominated = lambda a_list, b_list : not True in [a < b for a, b in zip(a_list, b_list)]
    is_equivalent = lambda a_list, b_list : not False in [a == b for a, b in zip(a_list, b_list)]
    pe_list = [True]*len(error_grid[0])
    for i in range(len(error_grid[0])):
        curr_error_list = [error_list[i] for error_list in error_grid]
        for j in range(len(error_grid[0])):
            other_error_list = [error_list[j] for error_list in error_grid]
            if not is_equivalent(curr_error_list, other_error_list) and is_dominated(curr_error_list, other_error_list):
                pe_list[i] = False
                break
    return pe_list

def find_knee_point(error_grid:list) -> int:
    """
    Identifies the knee point given a list of error lists

    Parameters:
    * `error_grid`: The list of list of errors

    Returns the index of the knee point
    """
    
    # Extract Pareto-efficient errors
    pe_list = find_pareto_efficiency(error_grid)
    pe_error_grid = [[error for error, pe in zip(error_list, pe_list) if pe] for error_list in error_grid]

    # Normalise the errors
    norm_error_grid = []
    for error_list in pe_error_grid:
        average_error = np.average(error_list)
        norm_error_list = [error/average_error for error in error_list]
        norm_error_grid.append(norm_error_list)
    
    # Comput distances to the ideal point (i.e., 0)
    distance_list = []
    for i in range(len(norm_error_grid[0])):
        square_list = [norm_error_list[i]**2 for norm_error_list in norm_error_grid]
        distance = math.sqrt(sum(square_list))
        distance_list.append(distance)
    
    # Return the index of the knee point
    min_distance = min(distance_list)
    pe_min_index = distance_list.index(min_distance)
    min_index = [i for i in range(len(pe_list)) if pe_list[i]][pe_min_index]
    return min_index

def get_termination(error_grid:list, thresholds:list):
    """
    Calculates the iteration that the algorithm should terminate

    Parameters:
    * `error_grid`: Lists of error lists

    Returns the termination iteration index (starting from 0)
    """
    
    # Identify indexes of errors below their thresholds
    bt_index_grid = []
    for error_list, threshold in zip(error_grid, thresholds):
        bt_index_list = [i for i in range(len(error_list)) if error_list[i] < threshold]
        bt_index_grid.append(bt_index_list)

    # Identify common indexes
    bt_set = set(bt_index_grid[0])
    for bt_index_list in bt_index_grid[1:]:
        bt_set.intersection_update(bt_index_list)
    
    # Identify earliest common index
    bt_list = sorted(list(bt_set))
    if bt_list == []:
        return len(error_grid[0])
    return bt_list[0]

def plot_errors(error_list:list, threshold:float) -> None:
    """
    Plots the errors

    Parameters:
    * `error_list`: List of errors
    * `threshold`:  Threshold value for errors
    """

    # Determine iterations
    iterations = [i+1 for i in range(len(error_list))]
    label_list = [2*i+1 for i in range(len(iterations)//2)]

    # Identify solutions below threshold
    bt_index_list = [i for i in range(len(error_list)) if error_list[i] < threshold]
    bt_iterations = [iterations[i] for i in bt_index_list]
    bt_error_list = [error_list[i] for i in bt_index_list]

    # Add general formatting
    plt.figure(figsize=(5,5), dpi=200)
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    
    # Add formatting for ticks nad labels
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Iterations", fontsize=14)
    plt.xlim(min(label_list)-1, max(label_list)+1)
    plt.xticks(ticks=label_list, labels=label_list)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)

    # Plot errors
    plt.plot(iterations, error_list, marker="o", linewidth=3, color=ERROR_COLOUR)
    plt.scatter(bt_iterations, bt_error_list, marker="o", s=8**2, edgecolor=BT_COLOUR, linewidth=2, color=ERROR_COLOUR, zorder=3)

    # Add legend
    handles = [
        plt.scatter([], [], marker="o", s=6**2, color=ERROR_COLOUR, label="Error"),
        plt.scatter([], [], marker="o", s=8**2, edgecolor=BT_COLOUR, linewidth=2, color="white", label="Below Threshold"),
    ]
    legend = plt.legend(handles=handles, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper right")
    plt.gca().add_artist(legend)

def get_errors_dict(sim_path:str):
    """
    Calculates the errors and stores it in a dictionary
    
    Parameters:
    * `sim_path`: Path to the simulation results
    """

    # Get experimental data
    exp_dict = csv_to_dict(EXP_DATA_PATH)
    max_strain = exp_dict["strain_intervals"][-1]
    eval_strains = list(np.linspace(0, max_strain, 50))

    # Get simulation results
    sim_dir_list  = [f"{sim_path}/{sim_dir}" for sim_dir in os.listdir(sim_path) if "simulate" in sim_dir]
    sum_path_list = [f"{sim_dir}/summary.csv" for sim_dir in sim_dir_list if os.path.exists(f"{sim_dir}/summary.csv")]
    sim_dict_list = [csv_to_dict(summary_path) for summary_path in sum_path_list]

    # Get parameters
    prm_dict_list = [read_params(f"{sim_dir}/params.txt") for sim_dir in sim_dir_list if os.path.exists(f"{sim_dir}/params.txt")]
    prm_vals_list = [[prm_dict[prm_name] for prm_name in prm_dict.keys() if prm_name in PARAM_NAMES] for prm_dict in prm_dict_list]

    # Get surrogate approximations
    sur_dir_list   = [f"{sim_path}/{dir_path}" for dir_path in os.listdir(sim_path) if "surrogate" in dir_path]
    sur_map_list   = [csv_to_dict(f"{sur_dir}/map.csv") for sur_dir in sur_dir_list]
    sur_model_list = [Model(f"{sur_dir}/sm.pt", f"{sur_dir}/map.csv", EXP_DATA_PATH, max_strain) for sur_dir in sur_dir_list]
    sur_dict_list  = [sur_model.get_response(prm_vals) for sur_model, prm_vals in zip(sur_model_list, prm_vals_list)]

    # Initialise error storage
    errors_dict = {"se": [], "ge": []}

    # Compare results
    for sim_dict, sur_dict, sur_map in zip(sim_dict_list, sur_dict_list, sur_map_list):

        # Identify grains the surrogate model can predict
        map_param_names = sur_map["param_name"] 
        grain_ids = [int(mpn.replace("g","").replace("_phi_1","")) for mpn in map_param_names if mpn.startswith("g") and mpn.endswith("_phi_1")]

        # Calculate stress error
        stress_error = get_stress(
            stress_list_1 = sur_dict["stress"],
            stress_list_2 = sim_dict["average_stress"],
            strain_list_1 = sur_dict["strain"],
            strain_list_2 = sim_dict["average_strain"],
            eval_strains  = eval_strains
        )

        # Calculate orientation error
        geodesic_grid = get_geodesics(
            grain_ids     = grain_ids,
            data_dict_1   = sim_dict,
            data_dict_2   = sur_dict,
            strain_list_1 = sim_dict["average_strain"],
            strain_list_2 = sur_dict["strain_intervals"],
            eval_strains  = eval_strains
        )
        geodesic_error = np.average([np.average(geodesic_list) for geodesic_list in geodesic_grid])

        # Add errors
        errors_dict["se"].append(stress_error)
        errors_dict["ge"].append(geodesic_error)

    # Return errors
    errors_dict["iteration"] = [int(i+1) for i in range(len(errors_dict["se"]))]
    return errors_dict

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
    sm_error()

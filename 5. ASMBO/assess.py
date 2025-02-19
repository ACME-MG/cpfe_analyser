"""
 Title:         Assesses the Pareto efficiency and tolerance
 Description:   Plots the errors for a set of simulations with Pareto efficiency and tolerance
 Author:        Janzen Choi

"""

# Libraries
import os
import matplotlib.pyplot as plt
import math, numpy as np
import sys; sys.path += [".."]
from __common__.io import csv_to_dict, dict_to_stdout
from __common__.general import round_sf
from __common__.analyse import get_geodesics, get_stress
from __common__.plotter import save_plot

# Paths
ASMBO_DIR     = "2025-02-14 (lh6_sm16_i34)" # 250213134340_i26_simulate
SIM_DATA_PATH = f"/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo/{ASMBO_DIR}"
EXP_DATA_PATH = "data/617_s3_40um_exp.csv"
RESULTS_PATH  = "results"

# Fields
STRAIN_FIELD = "average_strain"
STRESS_FIELD = "average_stress"

# Grain information
CAL_GRAIN_IDS = [51, 56, 72, 80, 126, 223, 237, 262]
VAL_GRAIN_IDS = [44, 50, 60, 178, 190, 207, 278, 299]

# Iteration information
END_ITER = 31
CUSTOM_LABEL_LIST = list(range(1,END_ITER+2,2)) # None

# Error information
ERROR_COLOUR = "black"
PD_COLOUR    = "tab:red"
PE_COLOUR    = "tab:green"
KP_COLOUR    = "tab:blue"
ACP_COLOUR   = "tab:green"
TOL_COLOUR   = "tab:blue"
SE_TOLERANCE = 0.02
GE_TOLERANCE = 0.08

def main():
    """
    Main function
    """
    tolerance()
    print()
    pareto()

def tolerance(sim_path:str=""):
    """
    Assesses a simulation based on its tolerance
    
    Parameters:
    * `sim_path`: Path to the simulation results
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
        plot_tolerance_error(se_list, SE_TOLERANCE)
        plt.ylabel(r"$E_{\sigma}$", fontsize=14)
        plt.yscale("log")
        # plt.ylim(10**-3, 10**9)
        save_plot("results/assess_tol_se.png")

        # Plot geodesic errors
        plot_tolerance_error(ge_list, GE_TOLERANCE)
        plt.ylabel(r"$E_{\phi}$", fontsize=14)
        plt.yscale("log")
        # plt.ylim(10**-3, 10**9)
        save_plot("results/assess_tol_ge.png")

        # Print out errors
        acceptability = [i in acceptable for i in range(len(se_list))]
        error_dict = {
            "iteration": list(range(1,len(se_list)+1)),
            "stress":    round_sf(se_list,5),
            "geodesic":  round_sf(ge_list,5),
            "accept":    acceptability,
        }
        print("Tolerance")
        dict_to_stdout(error_dict)

    # Otherwise, return the termination iteration
    else:
        if acceptable == []:
            return None
        return acceptable[0]+1

def plot_tolerance_error(error_list:list, tolerance:float) -> None:
    """
    Plots the errors with tolerance

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

    # Add legend to plot
    handles = [
        plt.scatter([], [], marker="o", s=6**2, color=ERROR_COLOUR, label="Error"),
        plt.scatter([], [], marker="o", s=8**2, color="white",      label="Under Tolerance",  edgecolor=TOL_COLOUR, linewidth=2),
    ]
    legend = plt.legend(handles=handles, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper right")
    plt.gca().add_artist(legend)

def pareto(sim_path:str=""):
    """
    Assesses a simulation based on its Pareto efficiency

    Parameters:
    * `sim_path`: Path to the simulation results
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

    # Calculate calibration error
    eval_strains = np.linspace(0, exp_dict["strain_intervals"][-1], 32)
    se_list, ge_list = get_errors(sim_dict_list, exp_dict, eval_strains, CAL_GRAIN_IDS)

    # Determine whether solutions are Pareto-efficient
    pe_list = find_pareto_efficiency([se_list, ge_list])

    # Determine whether solutions are Pareto-dominated
    pd_list = []
    for i in range(len(se_list)):
        dominated = False
        for prev_se, prev_ge in zip(se_list[:i], ge_list[:i]):
            if se_list[i] > prev_se and ge_list[i] > prev_ge:
                dominated = True
                break
        pd_list.append(dominated)

    # Determine the knee point
    kp_index = find_knee_point([se_list, ge_list])

    # If the initial simulation path was undefined, summarise results
    if sim_path == "":

        # Plot stress-strain error
        plot_pareto_error(se_list, pe_list, pd_list, kp_index)
        plt.ylabel(r"$E_{\sigma}$", fontsize=14)
        plt.yscale("log")
        # plt.ylim(10**-3, 10**9)
        save_plot("results/assess_par_se.png")

        # Plot geodesic error
        plot_pareto_error(ge_list, pe_list, pd_list, kp_index)
        plt.ylabel(r"$E_{\phi}$", fontsize=14)
        # plt.ylim(0, 0.7)
        save_plot("results/assess_par_ge.png")

        # Print out errors
        knee_point = [i==kp_index for i in range(len(se_list))]
        error_dict = {
            "iteration":  list(range(1,len(se_list)+1)),
            "stress":     round_sf(se_list,5),
            "geodesic":   round_sf(ge_list,5),
            "efficient":  pe_list,
            "dominated":  pd_list,
            "knee_point": knee_point,
        }
        print("Pareto Efficiency")
        dict_to_stdout(error_dict)

    # Otherwise, return the termination iteration
    else:
        termination = None
        for i in range(len(pd_list)-2):
            if pd_list[i] and pd_list[i+1] and pd_list[i+2]:
                termination = i+2+1 # +1 to start at 1
                break
        return termination

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

def plot_pareto_error(error_list:list, pe_list:list, pd_list:list, kp_index:int) -> None:
    """
    Plots the errors with Pareto-dominance

    Parameters:
    * `error_list`: The list of solution errors
    * `pe_list`:    The list of whether solutions are Pareto-efficient
    * `pd_list`:    The list of whether solutions are Pareto-dominated
    * `kp_index`:   The index of the knee point
    """

    # Plot the error
    label_list = list([i+1 for i in range(len(error_list))])
    initialise_error_plot(label_list)
    plt.plot(label_list, error_list, marker="o", linewidth=3, color=ERROR_COLOUR)

    # Highlight the Pareto-efficient errors
    pe_label_list = [label for label, pe in zip(label_list, pe_list) if pe]
    pe_error_list = [error for error, pe in zip(error_list, pe_list) if pe]
    plt.scatter(pe_label_list, pe_error_list, marker="o", s=8**2, edgecolor=PE_COLOUR, linewidth=2, color=ERROR_COLOUR, zorder=3)

    # Highlight the Pareto-dominated errors
    pd_label_list = [label for label, pd in zip(label_list, pd_list) if pd]
    pd_error_list = [error for error, pd in zip(error_list, pd_list) if pd]
    plt.scatter(pd_label_list, pd_error_list, marker="o", s=8**2, edgecolor=PD_COLOUR, linewidth=2, color=ERROR_COLOUR, zorder=3)

    # Highlight the knee point error
    plt.scatter([label_list[kp_index]], [error_list[kp_index]], marker="o", s=8**2, edgecolor=KP_COLOUR, linewidth=2, color=ERROR_COLOUR, zorder=3)

    # Add legend to plot
    handles = [
        plt.scatter([], [], marker="o", s=6**2, color=ERROR_COLOUR, label="Error"),
        plt.scatter([], [], marker="o", s=8**2, edgecolor=PD_COLOUR, linewidth=2, color="white", label="Previously dominated"),
        plt.scatter([], [], marker="o", s=8**2, edgecolor=PE_COLOUR, linewidth=2, color="white", label="Final Pareto front"),
        plt.scatter([], [], marker="o", s=8**2, edgecolor=KP_COLOUR, linewidth=2, color="white", label="Knee point"),
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
    plt.xlim(min(label_list)-1, max(label_list)+1)
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
    main()

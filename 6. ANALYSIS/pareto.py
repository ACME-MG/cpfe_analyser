"""
 Title:         Plot Optimised Simulations
 Description:   Plots the response of the optimised simulations
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += ["..", "/home/janzen/code/mms"]
import numpy as np
from __common__.io import csv_to_dict
from __common__.analyse import get_stress, get_geodesics
from __common__.general import round_sf, transpose
from __common__.pole_figure import get_lattice, IPF
from __common__.plotter import save_plot, prep_plot, set_limits, EXP_COLOUR
from __common__.orientation import rad_to_deg
import matplotlib.pyplot as plt
from __common__.surrogate import Model

# Colours
NORMAL_COLOUR = "darkgray"
EFFICIENT_COLOUR = "tab:green"
KNEEPOINT_COLOUR = "tab:red"

# Paths
EXP_DATA_PATH = "data/617_s3_40um_exp.csv"
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo/2025-03-18 (vh_x_sm8_i41)"
# PARAMS_PATH = f"{ASMBO_PATH}/250318013948_i21_optimise/params.csv"
PARAMS_PATH = "data/params.csv"
NUM_DATA = 100
# PARAMS_PATH = "data/params_all.csv"
# NUM_DATA = 200
SUR_PATH = f"{ASMBO_PATH}/250318012854_i21_surrogate"

def main():
    """
    Main function
    """

    # Get experimental data
    exp_dict = csv_to_dict(EXP_DATA_PATH)
    max_strain = exp_dict["strain_intervals"][-1]

    # Get parameters
    params_dict = csv_to_dict(PARAMS_PATH)
    param_fields = [param for param in params_dict.keys() if param.startswith("Param")]
    param_values = transpose([params_dict[pn][:NUM_DATA] for pn in param_fields])
    param_values += [[1938.04, 0.17367, 149.03, 6.1855]]

    # Evaluate surrogate model
    sur_model = Model(f"{SUR_PATH}/sm.pt", f"{SUR_PATH}/map.csv", EXP_DATA_PATH, max_strain)
    sur_dict_list = [sur_model.get_response(pv) for pv in param_values]

    # Identify grains the surrogate model can predict
    sur_map = csv_to_dict(f"{SUR_PATH}/map.csv")
    map_param_names = sur_map["param_name"] 
    grain_ids = [int(mpn.replace("g","").replace("_phi_1","")) for mpn in map_param_names if mpn.startswith("g") and mpn.endswith("_phi_1")]

    # Initialise errors
    eval_strains = list(np.linspace(0, max_strain, 50))
    stress_error_list = []
    geodesic_error_list = []

    # Obtain errors
    for sur_dict in sur_dict_list:
        stress_error, geodesic_error = get_errors(exp_dict, sur_dict, grain_ids, eval_strains)
        stress_error_list.append(stress_error)
        geodesic_error_list.append(geodesic_error)

    # Identify pareto front
    norm_se_list = normalise(stress_error_list)
    norm_ge_list = normalise(geodesic_error_list)
    error_grid = [norm_se_list, norm_ge_list]
    pe_list = find_pareto_efficiency(error_grid)

    # Identify pareto-efficient and non-efficient values
    ne_se_list = [norm_se_list[i] for i in range(len(pe_list)) if not pe_list[i]]
    ne_ge_list = [norm_ge_list[i] for i in range(len(pe_list)) if not pe_list[i]]
    pe_se_list = [norm_se_list[i] for i in range(len(pe_list)) if pe_list[i]]
    pe_ge_list = [norm_ge_list[i] for i in range(len(pe_list)) if pe_list[i]]
    distances = [pe_se**2 + pe_ge**2 for pe_se, pe_ge in zip(pe_se_list, pe_ge_list)]
    kp_index = distances.index(min(distances))

    # Adjust and plot errors
    prep_plot(r"$\hat{E}_{\phi}$", r"$\hat{E}_{\sigma}$", "", "")
    plt.scatter(ne_ge_list, ne_se_list, color=NORMAL_COLOUR, s=8**2)
    for i, (pe_se, pe_ge) in enumerate(zip(pe_se_list, pe_ge_list)):
        colour = EFFICIENT_COLOUR if i != kp_index else KNEEPOINT_COLOUR
        plt.plot([0, pe_ge], [0, pe_se], color=colour, linestyle="--")
        plt.scatter([pe_ge], [pe_se], color=colour, s=8**2)
    plt.scatter([0], [0], color="black", s=4**2, zorder=2)
    
    # Format error plot
    handles = [
        plt.scatter([], [], color=NORMAL_COLOUR,    label="Pareto-Dominated",  s=8**2),
        plt.scatter([], [], color=EFFICIENT_COLOUR, label="Pareto-Efficient",  s=8**2),
        plt.scatter([], [], color=KNEEPOINT_COLOUR, label="Knee Point",        s=8**2),
    ]
    legend = plt.legend(handles=handles, ncol=1, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    plt.gca().add_artist(legend)
    set_limits((-0.1, 1.1), (-0.1, 1.1))
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    save_plot("results/pareto_errors.png")

    # Get surrogate responses only for pareto efficient
    sur_dict_list = [sur_dict_list[i] for i in range(len(pe_list)) if pe_list[i]]

    # Initialise stress-strain plot
    prep_plot("Strain", "Stress", "mm/mm", "MPa")
    set_limits((0,0.5), (0,1800))

    # Plot stress-strain data
    plt.scatter(exp_dict["strain"], exp_dict["stress"], color=NORMAL_COLOUR, s=8**2)
    for i, sur_dict in enumerate(sur_dict_list):
        colour = EFFICIENT_COLOUR if i != kp_index else KNEEPOINT_COLOUR
        zorder = 4 if i != kp_index else 5
        plt.plot(sur_dict["strain"], sur_dict["stress"], color=colour, alpha=1.0, linewidth=3, zorder=zorder)

    # Format and save
    add_legend()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    save_plot("results/pareto_ss.png")

    # Initialise IPF
    ipf = IPF(get_lattice("fcc"))
    direction = [1,0,0]
    get_trajectories = lambda dict, g_ids : [transpose([dict[f"g{g_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]) for g_id in g_ids]

    # Plot experimental reorientation trajectories
    exp_trajectories = get_trajectories(exp_dict, grain_ids)
    if 14 in grain_ids:
        exp_trajectories[grain_ids.index(14)] = exp_trajectories[grain_ids.index(14)][:-1] 
    if 287 in grain_ids:
        exp_trajectories[grain_ids.index(287)] = exp_trajectories[grain_ids.index(287)][:-3] 
    ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": EXP_COLOUR, "linewidth": 3})
    ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": EXP_COLOUR, "head_width": 0.01, "head_length": 0.015})
    ipf.plot_ipf_trajectory([[et[0]] for et in exp_trajectories], direction, "scatter", {"color": EXP_COLOUR, "s": 8**2})

    # Iterate through simulations
    for i, sur_dict in enumerate(sur_dict_list):
        sim_trajectories = get_trajectories(sur_dict, grain_ids)
        if 44 in grain_ids:
            sim_trajectories[grain_ids.index(44)] = sim_trajectories[grain_ids.index(44)][:-9] 
        colour = EFFICIENT_COLOUR if i != kp_index else KNEEPOINT_COLOUR
        zorder = 4 if i != kp_index else 5
        ipf.plot_ipf_trajectory(sim_trajectories, direction, "plot", {"color": colour, "linewidth": 2, "zorder": zorder, "alpha": 1})
        ipf.plot_ipf_trajectory(sim_trajectories, direction, "arrow", {"color": colour, "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": zorder, "alpha": 1})
        ipf.plot_ipf_trajectory([[st[0]] for st in sim_trajectories], direction, "scatter", {"color": colour, "s": 6**2, "zorder": zorder, "alpha": 1})

    # Format and save IPF plot
    add_legend()
    save_plot(f"results/pareto_rt.png")

def add_legend():
    """
    Adds a legend to the current axis
    """
    handles = [
        plt.scatter([], [], color=EXP_COLOUR,    label="Experimental",     s=8**2),
        plt.plot([], [], color=EFFICIENT_COLOUR, label="Pareto-Efficient", linewidth=3)[0],
        plt.plot([], [], color=KNEEPOINT_COLOUR, label="Knee Point",       linewidth=3)[0],
    ]
    legend = plt.legend(handles=handles, ncol=1, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    plt.gca().add_artist(legend)

def normalise(value_list:list, min_norm:float=0.0, max_norm:list=1.0) -> list:
    """
    Normalises a list of values

    Parameters:
    * `value_list`: The list of values
    * `min_norm`:   The minimum value for the normalised list of values
    * `max_norm`:   The maximum value for the normalised list of values

    Returns the normalised list
    """
    min_value = min(value_list)
    max_value = max(value_list)
    normalised = [min_norm+((value-min_value)/(max_value-min_value))*(max_norm-min_norm) for value in value_list]
    return normalised

def get_errors(exp_dict:dict, sur_dict:dict, grain_ids:list, eval_strains:list) -> tuple:
    """
    Gets the errors

    Parameters
    * `exp_dict`:     Experimental data
    * `sur_dict`:     Surrogate model data
    * `grain_ids`:    List of grain IDs
    * `eval_strains`: Strains to evaluate errors

    Returns stress and geodesic errors as a tuple
    """

    # Calculate stress error
    stress_error = get_stress(
        stress_list_1 = sur_dict["stress"],
        stress_list_2 = exp_dict["stress"],
        strain_list_1 = sur_dict["strain"],
        strain_list_2 = exp_dict["strain"],
        eval_strains  = eval_strains
    )

    # Calculate orientation error
    geodesic_grid = get_geodesics(
        grain_ids     = grain_ids,
        data_dict_1   = exp_dict,
        data_dict_2   = sur_dict,
        strain_list_1 = exp_dict["strain_intervals"],
        strain_list_2 = sur_dict["strain_intervals"],
        eval_strains  = eval_strains
    )
    geodesic_error = np.average([np.average(geodesic_list) for geodesic_list in geodesic_grid])

    # Return errors
    return stress_error, geodesic_error

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

# Call main function
if __name__ == "__main__":
    main()

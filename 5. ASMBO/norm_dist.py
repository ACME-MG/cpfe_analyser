"""
 Title:         Normal Distributions
 Description:   Plots the errors of the optimised simulations using normal distributions
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats 
from __common__.io import csv_to_dict
from __common__.general import round_sf, pad_to_length
from __common__.plotter import save_plot, Plotter
from __common__.analyse import get_geodesics, get_stress

# Experimental Information
EXP_PATH = "data/617_s3_40um_exp.csv"
EXP_COLOUR = "silver"
EXP_EBSD_ID = "ebsd_4"

# Simulation Information
RES_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo"
SIM_INFO_LIST = [
    {"label": "VH",   "ebsd_id": "ebsd_4", "colour": "tab:cyan",   "path": "2025-01-18 (vh_sm8_i24)/250118195346_i8_simulate"},
    {"label": "LH-2", "ebsd_id": "ebsd_4", "colour": "tab:orange", "path": "2025-01-09 (lh_sm32_i16)/250108194247_i8_simulate"},
    {"label": "LH-6", "ebsd_id": "ebsd_4", "colour": "tab:purple", "path": "2025-01-18 (lh6_sm72_i20)/250117013234_i11_simulate"},
]
for sim_info in SIM_INFO_LIST:
    sim_info["data"] = csv_to_dict(f"{RES_PATH}/{sim_info['path']}/summary.csv")

# Other Constants
STRAIN_FIELD = "average_strain"
STRESS_FIELD = "average_stress"
RES_DATA_MAP = "data/res_grain_map.csv"

# Main function
def main():

    # Read experimental and simulated data
    exp_dict = csv_to_dict(EXP_PATH)
    eval_strains = exp_dict["strain_intervals"]
    sim_dict_list = [sim_info["data"] for sim_info in SIM_INFO_LIST]

    # Get the geodesic errors
    base_sim = SIM_INFO_LIST[0]["data"]
    grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in base_sim.keys() if "_phi_1" in key]
    ge_grid = get_geodesic_errors(sim_dict_list, exp_dict, eval_strains, grain_ids)

    # Initialise the plotter
    plotter = Plotter("strain", "", "mm/mm")
    plotter.prep_plot(size=14)
    plt.xlabel(r"$E_{\phi}$", fontsize=14)
    plt.ylabel("Probability Density", fontsize=14)
    plotter.set_limits((0,0.25), (0,20))
    plt.yticks([0, 5, 10, 15, 20])

    # Calculate normal distributions for each simulation
    for sim_info, ge_list in zip(SIM_INFO_LIST, ge_grid):

        # Get mean and standard deviation
        mean = np.average(ge_list)
        std  = np.std(ge_list)

        # Plot the distribution
        ge_x_list = np.linspace(min(ge_list)-1, max(ge_list)+1, 1000)
        ge_y_list = stats.norm.pdf(ge_x_list, mean, std)
        plt.plot(ge_x_list, ge_y_list, linewidth=3, color=sim_info["colour"])
        
        # Store distribution information
        mean_str = pad_to_length(round_sf(mean, 3), 6)
        std_str  = pad_to_length(round_sf(std, 3), 6)
        sim_info["norm"] = f"({mean_str}, {std_str})"
    
    # Define legend with supplementary information
    handles  = [plt.plot([], [], color=sim_info["colour"], label=sim_info["label"], marker="o", linewidth=3)[0] for sim_info in SIM_INFO_LIST]
    handles += [plt.scatter([], [], color="white", label=sim_info["norm"], marker="o", s=0) for sim_info in SIM_INFO_LIST]
    legend = plt.legend(handles=handles, ncol=2, columnspacing=-2.5, framealpha=1, edgecolor="black",
                        fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    plt.gca().add_artist(legend)

    # Format and save
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    save_plot("results/norm_dist_ge.png")

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
    plt.xticks(ticks=label_list, labels=label_list)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)

def get_stress_errors(sim_dict_list:list, exp_dict:dict, eval_strains:list) -> list:
    """
    Calculates the stress errors of a list of simulations relative to experimental data

    Parameters:
    * `sim_dict_list`: The list of dictionaries of simulation results
    * `exp_dict`:      The dictionary of experimental data
    * `eval_strains`:  The strains to conduct the error evaluations
    
    Returns the stress errors as a list
    """
    stress_error_list = []
    for sim_dict in sim_dict_list:
        stress_error = get_stress(
            stress_list_1 = exp_dict["stress"],
            stress_list_2 = sim_dict["average_stress"],
            strain_list_1 = exp_dict["strain"],
            strain_list_2 = sim_dict["average_strain"],
            eval_strains  = eval_strains
        )
        stress_error_list.append(stress_error)
    return stress_error_list

def get_geodesic_errors(sim_dict_list:list, exp_dict:dict, eval_strains:list, grain_ids:list) -> tuple:
    """
    Calculates the errors of a list of simulations relative to experimental data

    Parameters:
    * `sim_dict_list`: The list of dictionaries of simulation results
    * `exp_dict`:      The dictionary of experimental data
    * `eval_strains`:  The strains to conduct the error evaluations
    * `grain_ids`:     The list of grain IDs
    
    Returns the geodesic errors as a list
    """
    geodesic_errors_list = []
    for sim_dict in sim_dict_list:
        geodesic_grid = get_geodesics(
            grain_ids     = grain_ids,
            data_dict_1   = sim_dict,
            data_dict_2   = exp_dict,
            strain_list_1 = sim_dict["average_strain"],
            strain_list_2 = exp_dict["strain_intervals"],
            eval_strains  = eval_strains
        )
        geodesic_errors = [np.average(geodesic_list) for geodesic_list in geodesic_grid]
        geodesic_errors_list.append(geodesic_errors)
    return geodesic_errors_list

def get_sim_grain_id(exp_grain_id:int, ebsd_id:str) -> int:
    """
    Maps the experimental to simulated grain ID

    Parameters:
    * `exp_grain_id`: The grain ID from the experimental data
    * `ebsd_id`:      The origin of the EBSD map used to run the simulation

    Returns the corresponding grain ID in the simulated data
    """
    grain_map = csv_to_dict(RES_DATA_MAP)
    exp_ebsd_ids = grain_map[EXP_EBSD_ID]
    sim_ebsd_ids = grain_map[ebsd_id]
    return int(sim_ebsd_ids[exp_ebsd_ids.index(exp_grain_id)])

# Calls the main function
if __name__ == "__main__":
    main()

"""
 Title:         Normal Distributions
 Description:   Plots the errors of the optimised simulations using normal distributions
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats 
from __common__.io import csv_to_dict
from __common__.general import round_sf, pad_to_length
from __common__.plotter import save_plot, Plotter
from __common__.analyse import get_geodesics

# Experimental Information
EXP_PATH    = "data/617_s3_40um_exp.csv"
EXP_COLOUR  = "silver"
EXP_EBSD_ID = "ebsd_4"
# GRAIN_IDS   = [
#     51, 56, 126, 237, 262, 44, 60, 78, 86, 190, 53, 54, 59, 69, 82, 173, 254, 256, 271,
#     283, 303, 72, 80, 223, 178, 207, 244, 63, 77, 101, 117, 120, 141, 149, 157, 159, 193,
#     242, 255, 264, 273, 276, 277, 278, 280, 281, 286, 295, 299, 10, 37, 47, 76, 90, 91,
#     192, 202, 230, 235, 238, 259, 265, 306, 11, 20, 38, 39, 40, 64, 87, 107, 111, 128,
#     18, 85, 217, 284, 285
# ]
GRAIN_IDS = []
REMOVE_OUTLIERS = True

# Simulation Information
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo"
MOOSE_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim"
SIM_INFO_LIST = [
    # {"label": "Low-fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-02-02 (vh_sm8_i72)/250202092030_i59_simulate"},
    # {"label": "High-fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/2025-02-04 (617_s3_10um_vh)"},
    {"label": "Low-fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-02-06 (lh2_sm8_i15)/250206013451_i12_simulate"},
    {"label": "High-fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/2025-02-09 (617_s3_10um_lh2)"},
    # {"label": "VH",  "ebsd_id": "ebsd_2", "colour": "tab:cyan",   "path": f"{MOOSE_PATH}/2025-02-04 (617_s3_10um_vh)"},
    # {"label": "LH2", "ebsd_id": "ebsd_4", "colour": "tab:orange", "path": f"{ASMBO_PATH}/2025-01-09 (lh_sm32_i16)/250108194247_i8_simulate"},
    # {"label": "LH6", "ebsd_id": "ebsd_4", "colour": "tab:purple", "path": f"{ASMBO_PATH}/2025-01-18 (lh6_sm72_i20)/250117013234_i11_simulate"},
]
for si in SIM_INFO_LIST:
    si["data"] = csv_to_dict(f"{si['path']}/summary.csv")

# Other Constants
STRAIN_FIELD = "average_strain"
STRESS_FIELD = "average_stress"
RES_DATA_MAP = "data/res_grain_map.csv"
# SPACING      = -2.25
SPACING      = -2.25

# Main function
def main():

    # Read experimental and simulated data
    exp_dict = csv_to_dict(EXP_PATH)
    eval_strains = exp_dict["strain_intervals"]

    # Get the geodesic errors
    base_sim = SIM_INFO_LIST[0]["data"]
    grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in base_sim.keys() if "_phi_1" in key]
    grain_ids = grain_ids if GRAIN_IDS == [] else GRAIN_IDS
    ge_grid = get_geodesic_errors(SIM_INFO_LIST, exp_dict, eval_strains, grain_ids)

    # Initialise the plotter
    plotter = Plotter("strain", "", "mm/mm")
    plotter.prep_plot(size=14)
    plt.xlabel(r"$E_{\phi}$", fontsize=14)
    plt.ylabel("Probability Density", fontsize=14)
    plotter.set_limits((0,0.25), (0,25))
    plt.yticks([0, 5, 10, 15, 20, 25])

    # Calculate normal distributions for each simulation
    for sim_info, ge_list in zip(SIM_INFO_LIST, ge_grid):

        # Get mean and standard deviation
        if REMOVE_OUTLIERS:
            ge_list = remove_outliers(ge_list)
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
    legend = plt.legend(handles=handles, ncol=2, columnspacing=SPACING, framealpha=1, edgecolor="black",
                        fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    plt.gca().add_artist(legend)

    # Format and save
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    save_plot("results/norm_dist_ge.png")

def remove_outliers(data_list:list):
    """
    Removes outliers from a dataset using the IQR method
    
    Parameters:
    * `data_list`: List of data values

    Returns a new list without outliers
    """
    data_series = pd.Series(data_list)
    q1 = data_series.quantile(0.25)
    q3 = data_series.quantile(0.75)
    lower_bound = q1 - 1.5*(q3-q1)
    upper_bound = q3 + 1.5*(q3-q1)
    new_list = data_series[(data_series >= lower_bound) & (data_series <= upper_bound)].tolist()
    return new_list

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

def get_geodesic_errors(sim_info_list:list, exp_dict:dict, eval_strains:list, grain_ids:list) -> tuple:
    """
    Calculates the errors of a list of simulations relative to experimental data

    Parameters:
    * `sim_info_list`: The list of dictionaries of simulation results
    * `exp_dict`:      The dictionary of experimental data
    * `eval_strains`:  The strains to conduct the error evaluations
    * `grain_ids`:     The list of grain IDs
    
    Returns the geodesic errors as a list
    """

    # Iterate through simulations
    geodesic_error_list = []
    for si in sim_info_list:

        # Convert grain IDs
        eval_grain_ids = []
        sim_grain_ids = [get_sim_grain_id(grain_id, si["ebsd_id"]) for grain_id in grain_ids]
        sim_dict = {}
        for grain_id, sim_grain_id in zip(grain_ids, sim_grain_ids):
            if sim_grain_id == -1:
                continue
            for phi in ["phi_1", "Phi", "phi_2"]:
                sim_dict[f"g{grain_id}_{phi}"] = si["data"][f"g{sim_grain_id}_{phi}"]
            eval_grain_ids.append(grain_id)

        # Calculate geodesic errors
        geodesic_grid = get_geodesics(
            grain_ids     = eval_grain_ids,
            data_dict_1   = sim_dict,
            data_dict_2   = exp_dict,
            strain_list_1 = si["data"]["average_strain"],
            strain_list_2 = exp_dict["strain_intervals"],
            eval_strains  = eval_strains
        )

        # Compile geodesic errors
        geodesic_errors = [np.average(geodesic_list) for geodesic_list in geodesic_grid]
        geodesic_error_list.append(geodesic_errors)
    
    # Return
    return geodesic_error_list

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

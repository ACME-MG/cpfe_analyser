"""
 Title:         Plot Strain
 Description:   Plots the errors of the optimised simulations based on each strain
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from __common__.io import csv_to_dict
from __common__.general import transpose
from __common__.plotter import save_plot, Plotter
from __common__.analyse import get_geodesics

# Experimental Information
EXP_PATH    = "data/617_s3_40um_exp.csv"
EXP_COLOUR  = "silver"
EXP_EBSD_ID = "ebsd_4"
# GRAIN_IDS   = [51, 56, 72, 80, 126, 223, 237, 262]  # Calibration
# GRAIN_IDS   = [44, 50, 60, 178, 190, 207, 278, 299] # Validation
GRAIN_IDS   = []
REMOVE_OUTLIERS = True

# Simulation Information
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo"
MOOSE_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim"
BASE_SIM   = csv_to_dict(f"{ASMBO_PATH}/2025-02-02 (vh_sm8_i72)/250202092030_i59_simulate/summary.csv")
SIM_INFO_LIST = [

    # Voce Hardening Model
    # {"label": "Low-fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-02-02 (vh_sm8_i72)/250202092030_i59_simulate"},
    # {"label": "High-fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/2025-02-04 (617_s3_10um_vh)"},
    
    # Two-coefficient latent hardening model
    # {"label": "Low-fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-02-06 (lh2_sm8_i15)/250206013451_i12_simulate"},
    # {"label": "High-fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/2025-02-09 (617_s3_10um_lh2)"},
    
    # Six-coefficient latent hardening model
    # {"label": "Low-fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-02-14 (lh6_sm16_i34)/250213134340_i26_simulate"},
    # {"label": "High-fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/2025-02-14 (617_s3_10um_lh6)"},
    
    # All low-fidelity
    # {"label": "VH",   "ebsd_id": "ebsd_4", "colour": "tab:cyan",   "path": f"{ASMBO_PATH}/2025-02-02 (vh_sm8_i72)/250202092030_i59_simulate"},
    # {"label": "LH2",  "ebsd_id": "ebsd_4", "colour": "tab:orange", "path": f"{ASMBO_PATH}/2025-02-06 (lh2_sm8_i15)/250206013451_i12_simulate"},
    # {"label": "LH6",  "ebsd_id": "ebsd_4", "colour": "tab:purple", "path": f"{ASMBO_PATH}/2025-02-14 (lh6_sm16_i34)/250213134340_i26_simulate"},
    
    # All high-fidelity
    {"label": "VH",   "ebsd_id": "ebsd_2", "colour": "tab:cyan",   "path": f"{MOOSE_PATH}/2025-02-04 (617_s3_10um_vh)"},
    {"label": "LH2",  "ebsd_id": "ebsd_2", "colour": "tab:orange", "path": f"{MOOSE_PATH}/2025-02-09 (617_s3_10um_lh2)"},
    {"label": "LH6",  "ebsd_id": "ebsd_2", "colour": "tab:purple", "path": f"{MOOSE_PATH}/2025-02-14 (617_s3_10um_lh6)"},
]
for si in SIM_INFO_LIST:
    si["data"] = csv_to_dict(f"{si['path']}/summary.csv")

# Other Constants
STRAIN_FIELD = "average_strain"
STRESS_FIELD = "average_stress"
RES_DATA_MAP = "data/res_grain_map.csv"

# Main function
def main():

    # Read experimental and simulated data
    exp_dict = csv_to_dict(EXP_PATH)
    eval_strains = exp_dict["strain_intervals"]

    # Get the geodesic errors
    grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in BASE_SIM.keys() if "_phi_1" in key]
    grain_ids = grain_ids if GRAIN_IDS == [] else GRAIN_IDS

    ge_grid = get_geodesic_errors(SIM_INFO_LIST, exp_dict, eval_strains, grain_ids)

    # Initialise the plotter
    plotter = Plotter("strain", "", "mm/mm")
    plotter.prep_plot(size=14)
    plt.xlabel("Strain (mm/mm)", fontsize=14)
    plt.ylabel(r"$E_{\phi}$", fontsize=14)
    # plotter.set_limits((0,0.25), (0,25))
    # plt.yticks([0, 5, 10, 15, 20, 25])

    # Plot each simulation result
    for sim_info, ge_list in zip(SIM_INFO_LIST, ge_grid):
        strain_list = eval_strains
        if REMOVE_OUTLIERS:
            new_ge_list = remove_outliers(ge_list)
            strain_list = [strain for strain, ge in zip(strain_list, ge_list) if ge in new_ge_list]
            ge_list = new_ge_list
        plt.plot(strain_list, ge_list, color=sim_info["colour"], label=sim_info["label"], marker="o", linewidth=3)
    
    # Define legend with supplementary information
    handles  = [plt.plot([], [], color=sim_info["colour"], label=sim_info["label"], marker="o", linewidth=3)[0] for sim_info in SIM_INFO_LIST]
    legend = plt.legend(handles=handles, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    plt.gca().add_artist(legend)

    # Format and save
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    save_plot("results/strain_ge.png")

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
    Calculates the errors of a list of simulations relative to experimental data for each strain

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
        geodesic_errors = [np.average(geodesic_list) for geodesic_list in transpose(geodesic_grid)]
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
    if exp_grain_id in exp_ebsd_ids:
        return int(sim_ebsd_ids[exp_ebsd_ids.index(exp_grain_id)])
    else:
        return None

# Calls the main function
if __name__ == "__main__":
    main()

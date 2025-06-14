"""
 Title:         Plot Reorientation Trajectory Errors
 Description:   Plots the errors of individual reorientation trajectories
 Author:        Janzen Choi

"""

# Libraries
import os
import matplotlib.pyplot as plt
import math, numpy as np
import sys; sys.path += [".."]
from __common__.analyse import get_geodesics
from __common__.io import csv_to_dict
from __common__.plotter import save_plot

# Experimental Information
EXP_DATA_PATH = "data/617_s3_40um_exp.csv"
EXP_EBSD_ID   = "ebsd_4"
RES_DATA_MAP  = "data/res_grain_map.csv"

# Simulation Information
RESULTS_PATH = "results"
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo"
MOOSE_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/moose_sim"
SIM_DATA_LIST = [
    {"label": "VH",  "path": f"{MOOSE_PATH}/2025-03-15 (617_s3_vh_x_hr)",     "ebsd_id": "ebsd_2", "colour": "tab:cyan"},
    {"label": "LH2", "path": f"{MOOSE_PATH}/2025-04-05 (617_s3_lh2_di_x_hr)", "ebsd_id": "ebsd_2", "colour": "tab:orange"},
    {"label": "LH6", "path": f"{MOOSE_PATH}/2025-04-28 (617_s3_lh6_di_x_hr)", "ebsd_id": "ebsd_2", "colour": "tab:purple"},
]

# Grain IDs
GRAIN_IDS = [
    [14, 72, 95, 101, 207, 240, 262, 287],  # Calibration
    [39, 50, 138, 164, 185, 223, 243, 238], # Validation
]

# Plotting Information
MAX_ERROR = 0.30

def main():
    """
    Main function
    """
    
    # Iterate through grains
    offset = 0
    for i, grain_ids in enumerate(GRAIN_IDS):
        
        # Initialise plot
        prep_plot("Grain IDs", r"$E_{\Phi}$")
        
        # Offset grain IDs
        new_grain_ids = [i+1+offset for i in range(len(grain_ids))]
        offset += len(grain_ids)
        set_limits((min(new_grain_ids)-0.5, max(new_grain_ids)+0.5), (0, MAX_ERROR))

        # Plot errors for each simulation
        for sim_data in SIM_DATA_LIST:
            error_list = get_error_list(sim_data["path"], sim_data["ebsd_id"], grain_ids)
            print(f"Min: {min(error_list)}\tMax: {max(error_list)}")
            plt.scatter(new_grain_ids, error_list, marker="o", s=8**2, edgecolor="black", linewidth=1, color=sim_data["colour"], zorder=3)
        
        # Add legend
        label_list = [sim_data["label"] for sim_data in SIM_DATA_LIST]
        colour_list = [sim_data["colour"] for sim_data in SIM_DATA_LIST]
        add_legend(label_list, colour_list)

        # Save plot
        save_plot(f"{RESULTS_PATH}/plot_rt_error_{i+1}.png")

def prep_plot(x_label:str, y_label:str, title:str="", size:int=14) -> None:
    """
    Prepares the plot
    
    Parameters:
    * `x_label`: Label for the x-axis
    * `y_label`: Label for the y-axis
    * `title`:   The title of the plot
    * `size`:    The size of the font
    """

    # Set figure size and title
    plt.figure(figsize=(5,5), dpi=200)
    plt.title(title, fontsize=size+3, fontweight="bold", y=1.05)
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)

    # Set x and y labels
    plt.xlabel(x_label, fontsize=size)
    plt.ylabel(y_label, fontsize=size)
    
    # Format and save
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)

def set_limits(x_limits:tuple=None, y_limits:tuple=None) -> None:
    """
    Sets the limits of the x and y scales

    Parameters:
    * `x_limits`: The upper and lower bounds of the plot for the x scale
    * `y_limits`: The upper and lower bounds bound of the plot for the y scale
    """
    if x_limits != None:
        plt.xlim(*x_limits)
    if y_limits != None:
        plt.ylim(*y_limits)

def add_legend(label_list:list, colour_list:list) -> None:
    """
    Adds a basic legend

    Parameters:
    * `label_list`:  List of labels
    * `colour_list`: List of colours
    """
    handles = [plt.scatter([], [], color=colour, label=label, marker="o", s=8**2, edgecolor="black", linewidth=1)
               for label, colour in zip(label_list, colour_list)]
    legend = plt.legend(handles=handles, ncol=1, framealpha=1, edgecolor="black",
                        fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    plt.gca().add_artist(legend)

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

def get_error_list(sim_path:str, ebsd_id:str, grain_ids:list) -> list:
    """
    Calculates the errors in the reorientation trajectories
    
    Parameters:
    * `sim_path`:  Path to the simulation results
    * `ebsd_id`:   ID of the EBSD downsample
    * `grain_ids`: List of grain IDs to calculate the error

    Returns the list of errors
    """

    # Get experimental data
    exp_dict = csv_to_dict(EXP_DATA_PATH)
    max_strain = exp_dict["strain_intervals"][-1]
    eval_strains = list(np.linspace(0, max_strain, 50))
    
    # Get simulated data
    sim_dict = csv_to_dict(f"{sim_path}/summary.csv")
    sim_grain_ids = [get_sim_grain_id(grain_id, ebsd_id) for grain_id in grain_ids]

    # Convert the grain IDs
    new_sim_dict = {}
    for grain_id, sim_grain_id in zip(grain_ids, sim_grain_ids):
        for phi in ["phi_1", "Phi", "phi_2"]:
            new_sim_dict[f"g{grain_id}_{phi}"] = sim_dict[f"g{sim_grain_id}_{phi}"]

    # Calculate orientation error
    geodesic_grid = get_geodesics(
        grain_ids     = grain_ids,
        data_dict_1   = new_sim_dict,
        data_dict_2   = exp_dict,
        strain_list_1 = sim_dict["average_strain"],
        strain_list_2 = exp_dict["strain_intervals"],
        eval_strains  = eval_strains
    )

    # Format and return
    geodesic_error = [np.average(geodesic_list) for geodesic_list in geodesic_grid]
    return geodesic_error

# Calls the main function
if __name__ == "__main__":
    main()


"""
 Title:         Plot Texture
 Description:   Plots the texture of the optimised simulations at certain strains
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import numpy as np
import matplotlib.pyplot as plt
from __common__.general import transpose
from __common__.io import csv_to_dict, dict_to_csv
from __common__.plotter import save_plot, lighten_colour
from __common__.pole_figure import get_lattice, PFD
from __common__.interpolator import intervaluate

# Paths
EXP_PATH     = "data/617_s3_40um_exp.csv"
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo"
MOOSE_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/moose_sim"
RESULTS_PATH = "results/pf"

# Simulation information
SIM_INFO_LIST = [
    
    # VH model
    # {"label": "Low-Fidelity",  "colour": "tab:green", "lightness": (0.95, 0), "path": f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)/250318014435_i21_simulate"},
    # {"label": "High-Fidelity", "colour": "tab:red",   "lightness": (0.95, 0), "path": f"{MOOSE_PATH}/2025-03-15 (617_s3_vh_x_hr)"},
    
    # VH model alt
    # {"label": "Low-Fidelity",  "colour": "tab:green", "lightness": (0.95, 0), "path": f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i52_val)/250609222901_i8_simulate"},
    # {"label": "High-Fidelity", "colour": "tab:red",   "lightness": (0.95, 0), "path": f"{MOOSE_PATH}/2025-06-12 (617_s3_vh_di_x_hr_val)"},

    # LH2 model
    # {"label": "Low-Fidelity",  "colour": "tab:green", "lightness": (0.95, 0), "path": f"{ASMBO_PATH}/2025-03-28 (lh2_x_sm8_i29)/250327093649_i16_simulate"},
    # {"label": "High-Fidelity", "colour": "tab:red",   "lightness": (0.95, 0), "path": f"{MOOSE_PATH}/2025-04-05 (617_s3_lh2_di_x_hr)"},
    
    # LH6 model
    # {"label": "Low-Fidelity",  "colour": "tab:green", "lightness": (0.95, 0), "path": f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250422034348_i36_simulate"},
    # {"label": "High-Fidelity", "colour": "tab:red",   "lightness": (0.95, 0), "path": f"{MOOSE_PATH}/2025-04-28 (617_s3_lh6_di_x_hr)"},

    # All models (high-fidelity)
    # {"label": "VH",  "colour": "tab:cyan",   "lightness": (0.95, 0), "path": f"{MOOSE_PATH}/2025-03-15 (617_s3_vh_x_hr)"},
    # {"label": "LH2", "colour": "tab:orange", "lightness": (0.95, 0), "path": f"{MOOSE_PATH}/2025-04-05 (617_s3_lh2_di_x_hr)"},
    # {"label": "LH6", "colour": "tab:red",    "lightness": (0.95, 0), "path": f"{MOOSE_PATH}/2025-04-28 (617_s3_lh6_di_x_hr)"},
    
    # All models (low-fidelity)
    {"label": "VH",  "colour": "tab:cyan",   "lightness": (0.95, 0), "path": f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)/250318014435_i21_simulate"},
    {"label": "LH2", "colour": "tab:orange", "lightness": (0.95, 0), "path": f"{ASMBO_PATH}/2025-03-28 (lh2_x_sm8_i29)/250327093649_i16_simulate"},
    {"label": "LH6", "colour": "tab:red",    "lightness": (0.95, 0), "path": f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250422034348_i36_simulate"},
]
EXP_INFO = {"colour": "silver", "lightness": (0.95, 0)}

# Plotting parameters
DIRECTIONS   = [[1,1,1]]
# EVAL_STRAINS = [0.01, 0.052, 0.106, 0.208, 0.29]
EVAL_STRAINS = [0.0001, 0.00063414, 0.00153, 0.00494, 0.0098, 0.01483, 0.02085, 0.02646, 0.03516, 0.04409, 0.05197, 0.06013, 0.07059, 0.08208, 0.09406, 0.10561, 0.11929, 0.13656, 0.15442, 0.18237, 0.20849, 0.23627, 0.26264, 0.28965]
NUM_LEVELS   = 5
MRD_VALUES   = [0.0, 0.15, 0.30, 0.45, 0.60, 0.75] # broken

# Output parameters
SAVE_PLOT   = False
SAVE_LEGEND = False
SAVE_FILES  = True

# Formatting parameters
FONTSIZE     = 14
LINEWIDTH    = 1.5
FIGURE_X     = 5.2
FIGURE_Y     = 8
SCALE_FACTOR = 10/FIGURE_Y
DIM_RATIO    = (FIGURE_X/FIGURE_Y)
FIELD_LIST   = ["phi_1", "Phi", "phi_2", "volume", "area"]

def main():
    """
    Main function
    """

    # Initialise
    get_grain_ids    = lambda dict : [int(key.replace("g","").replace("_phi_1","")) for key in dict.keys() if "_phi_1" in key]
    get_trajectories = lambda dict, grain_ids : [transpose([dict[f"g{grain_id}_{phi}"] for phi in FIELD_LIST if f"g{grain_id}_{phi}" in dict.keys()]) for grain_id in grain_ids]

    # Read experimental data
    exp_dict = csv_to_dict(EXP_PATH)
    exp_strain_list = exp_dict["strain_intervals"]
    exp_grain_ids = get_grain_ids(exp_dict)
    exp_trajectories = get_trajectories(exp_dict, exp_grain_ids)

    # Read simulation data
    for sim_info in SIM_INFO_LIST:
        sim_info["dict"] = csv_to_dict(f"{sim_info['path']}/summary.csv")
        sim_info["strain"] = sim_info["dict"]["average_strain"]
        sim_grain_ids = get_grain_ids(sim_info["dict"])
        sim_info["trajectories"] = get_trajectories(sim_info["dict"], sim_grain_ids)

    # Extract orientation data
    exp_orientations_list = []
    sim_orientations_list_list = []
    for i, eval_strain in enumerate(EVAL_STRAINS):
        exp_orientations = [intervaluate_orientation(exp_strain_list, et, eval_strain) for et in exp_trajectories]
        sim_orientations_list = []
        for sim_info in SIM_INFO_LIST:
            sim_orientations = [intervaluate_orientation(sim_info["strain"], st, eval_strain) for st in sim_info["trajectories"]]
            sim_orientations_list.append(sim_orientations)
        exp_orientations_list.append(exp_orientations)
        sim_orientations_list_list.append(sim_orientations_list)

    # Generates plots
    if SAVE_PLOT:
        
        # Formatting settings
        figure_size = (5*DIM_RATIO, 5)
        top = 0.85
        right = 0.925
        left = 0.0
        adjust = {"top": top, "bottom": top-DIM_RATIO*(right-left), "right": right, "left": left}

        # Plot orientations for each evaluation strain
        for i, eval_strain in enumerate(EVAL_STRAINS):
            if i == 0:
                eval_strain = 0
            
            # Retrieve orientations
            exp_orientations = exp_orientations_list[i]
            sim_orientations_list = sim_orientations_list_list[i]

            # Plot orientations for each direction
            for direction in DIRECTIONS:

                # Initialise density-contoured pole figure
                pfd = PFD(get_lattice("fcc"), figsize=figure_size, adjust=adjust, fontsize=FONTSIZE*SCALE_FACTOR, linewidth=LINEWIDTH*SCALE_FACTOR)
                plt.text(0.0, -1.2, f"({round(eval_strain*100, 1)}%)", fontsize=FONTSIZE*SCALE_FACTOR, ha="center", va="center")

                # Plot experimental contours
                exp_lightness = np.linspace(EXP_INFO["lightness"][0], EXP_INFO["lightness"][1], NUM_LEVELS)
                exp_colours = [lighten_colour(EXP_INFO["colour"], el) for el in exp_lightness]
                pfd.plot_pfd(exp_orientations, direction, levels=NUM_LEVELS, colour_list=exp_colours)
                
                # Plot simulated contours
                for sim_orientations, sim_info in zip(sim_orientations_list, SIM_INFO_LIST):
                    sim_lightness = np.linspace(sim_info["lightness"][0], sim_info["lightness"][1], NUM_LEVELS)
                    sim_colours = [lighten_colour(sim_info["colour"], sl) for sl in sim_lightness]
                    pfd.plot_pfd(sim_orientations, direction, levels=NUM_LEVELS, colour_list=sim_colours)
                
                # Save the pole figure
                dir_str = "".join([str(d) for d in direction])
                save_plot(f"{RESULTS_PATH}/con_{dir_str}_e{i+1}")
            
    # Creates a legend
    if SAVE_LEGEND:
        colour_grid = [[lighten_colour(si["colour"], ln) for ln in np.linspace(si["lightness"][0], si["lightness"][1], len(MRD_VALUES))] for si in SIM_INFO_LIST]
        colour_grid = [[lighten_colour(EXP_INFO["colour"], ln) for ln in np.linspace(EXP_INFO["lightness"][0], EXP_INFO["lightness"][1], len(MRD_VALUES))]] + colour_grid
        add_legend(
            mrd_labels      = ["Experimental"] + [si["label"] for si in SIM_INFO_LIST],
            mrd_values      = MRD_VALUES,
            mrd_colour_grid = colour_grid,
            figure_factor   = 5/8,
        )
    
    # Exports orientation data
    if SAVE_FILES:
        for i in range(len(EVAL_STRAINS)):
            exp_orientations = exp_orientations_list[i]
            sim_orientations_list = sim_orientations_list_list[i]
            export_data(exp_orientations, f"exp_i{i+1}")
            for j, sim_orientations in enumerate(sim_orientations_list):
                export_data(sim_orientations, f"sim_s{j+1}_i{i+1}")

def export_data(euler_list:list, file_name:str) -> None:
    """
    Exports texture data as CSV files

    Parameters:
    * `euler_list`: The list of orientations in euler-bunge form (rads)
    * `file_name`:  Name of file to save the data (without extension)
    """
    phi_lists = transpose(euler_list)
    phi_dict = {}
    for i, field in enumerate(FIELD_LIST):
        if i < len(phi_lists):
            phi_dict[field] = phi_lists[i]
    dict_to_csv(phi_dict, f"{RESULTS_PATH}/{file_name}.csv")

def add_legend(mrd_labels:list, mrd_values:list, mrd_colour_grid:list, figure_factor:float=5/8) -> None:
    """
    Adds a legend containing the MRD key

    Parameters:
    * `mrd_labels`:      Labels for each model
    * `mrd_values`:      Levels of MRD
    * `mrd_colour_grid`: Grid of colours for the MRD values
    * `figure_factor`:   Factor to scale the figure's dimensions
    """

    # Initialise plot
    figure = plt.figure(figsize=(5*figure_factor,5), dpi=200)
    figure.subplots_adjust(top=0.83, bottom=0.266, left=0.34, right=0.85)

    # Draw bars
    for i, mrd_colour_list in enumerate(mrd_colour_grid):
        for j in range(len(mrd_values)):
            mrd_colour = mrd_colour_list[j]
            plt.bar(i+1, 1, color=mrd_colour, width=1, edgecolor="black", zorder=5, bottom=j)

    # Format grid
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(2*10/7)

    # Format ticks
    x_ticks = [i+1.3 for i in range(len(mrd_labels))]
    y_ticks = [i+0.5 for i in range(len(mrd_values))]
    plt.xticks(ticks=x_ticks, labels=mrd_labels, rotation=45, ha="right", fontsize=12*SCALE_FACTOR)
    plt.yticks(ticks=y_ticks, labels=mrd_values, fontsize=12*SCALE_FACTOR)
    plt.tick_params(axis="x", which="both", bottom=False, top=False)
    plt.tick_params(axis="y", which="both", left=False, right=False)

    # Apply labels and limits
    plt.ylabel("MRD", fontsize=14*SCALE_FACTOR, labelpad=14)
    plt.xlim(0.5, len(mrd_colour_grid)+0.5)
    plt.ylim(0, len(mrd_values))

    # Save
    save_plot(f"{RESULTS_PATH}/legend.png")

def intervaluate_orientation(strain_list:list, orientation_list:list, strain:float) -> list:
    """
    Interpolates and evaluates from a list of orientations

    Parameters:
    * `strain_list`:      List of strain values
    * `orientation_list`: List of orientation components
    * `strain`:           Strain value to evaluate

    Returns the intervaluated orientation
    """
    components = transpose(orientation_list)
    new_components = [intervaluate(strain_list, component, strain) for component in components]
    return transpose(new_components)

# Main function
if __name__ == "__main__":
    main()

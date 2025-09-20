"""
 Title:         Plot Auxiliary Texture Stuff
 Description:   Plots the texture of the optimised simulations at certain strains
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import numpy as np
import matplotlib.pyplot as plt
from __common__.io import csv_to_dict
from __common__.plotter import save_plot, Plotter

# Colours
EXP_COLOUR = {"base": "silver",     "start": [.95, .95, .95], "end": [0.6, 0.6, 0.6]}
VAL_COLOUR = {"base": "tab:red",    "start": [1.0, 0.9, 0.9], "end": [0.7, 0.1, 0.1]}
CAL_COLOUR = {"base": "tab:green",  "start": [.85, 0.9, .85], "end": [0.0, 0.5, 0.1]}
VH_COLOUR  = {"base": "tab:cyan",   "start": [0.9, 1.0, 1.0], "end": [0.0, 0.8, 0.8]}
LH2_COLOUR = {"base": "tab:orange", "start": [1.0, 0.9, 0.8], "end": [1.0, 0.5, 0.2]}
LH6_COLOUR = {"base": "tab:purple", "start": [0.9, 0.8, 1.0], "end": [0.5, 0.0, 0.8]}
OLV_COLOUR = {"base": "tab:olive",  "start": [1.0, 1.0, 0.9], "end": [0.8, 0.8, 0.0]}
BLU_COLOUR = {"base": "tab:blue",   "start": [0.9, 0.9, 1.0], "end": [0.1, 0.1, 0.7]}

# Simulation parameters
# HEADERS, KP_INDEX, COLOURS = ["Experimental"] + [f"Run {i+1}" for i in range(5)], 3, [EXP_COLOUR]+5*[CAL_COLOUR] # Runs
# COLOURS[KP_INDEX+1] = VAL_COLOUR
# HEADERS, KP_INDEX, COLOURS = ["Experimental", "Low-Fidelity", "High-Fidelity"], 1, [EXP_COLOUR, VAL_COLOUR, BLU_COLOUR] # Fidelities
HEADERS, KP_INDEX, COLOURS = ["Experimental", "VH", "LH2", "LH6"], 1, [EXP_COLOUR, VH_COLOUR, LH2_COLOUR, LH6_COLOUR] # Models

# Paths
DATA_PATH    = "/mnt/c/Users/janzen/Desktop/texture/scripts/results/texture_indexes.csv"
# DATA_PATH    = "data/texture_indexes.csv"
RESULTS_PATH = "results"

# Formatting parameters
FONTSIZE     = 14
TI_SCALE     = 5/4
FIGURE_X     = 5*(3/10)
FIGURE_Y     = 5*(8/10)
# FIGURE_X     = 5*(10/10)
# FIGURE_Y     = 5*(2.5/10)
SCALE_FACTOR = 10/FIGURE_Y

# Alpha parameters (None to deactivate)
ALPHA_LIST = None

# Plotting constants
EVAL_STRAINS = [0, 0.00063414, 0.00153, 0.00494, 0.0098, 0.01483, 0.02085, 0.02646, 0.03516, 0.04409, 0.05197, 0.06013, 0.07059, 0.08208, 0.09406, 0.10561, 0.11929, 0.13656, 0.15442, 0.18237, 0.20849, 0.23627, 0.26264, 0.28965]
LEVELS       = [0, 1, 2, 3, 4, 5]

def main():
    """
    Main function
    """
    
    # Plot texture indexes
    plot_ti(
        ti_dict    = csv_to_dict(DATA_PATH),
        ti_colours = [c["base"] for c in COLOURS]
    )
    
    # # Plot MRD colour chart
    # colour_grid = [np.linspace(c["start"], c["end"], len(LEVELS)) for c in COLOURS]
    # plot_mrd(
    #     mrd_labels      = ["Exp.", "LF ", "HF "],
    #     mrd_values      = LEVELS,
    #     mrd_colour_grid = colour_grid
    # )

def plot_ti(ti_dict:list, ti_colours:list) -> None:
    """
    Plots the texture indexes

    Parameters:
    * `ti_dict`:    Dictionary of texture indexes
    * `ti_colours`: List of colours
    """

    # Prepare texture indexes
    ti_grid = [ti_dict[header] for header in HEADERS]

    # Initialise plot
    plotter = Plotter()
    plotter.prep_plot(size=FONTSIZE*TI_SCALE)
    plt.xlabel("Strain (mm/mm)", fontsize=FONTSIZE*TI_SCALE)
    plt.ylabel("Texture Indexes", fontsize=FONTSIZE*TI_SCALE)
    plotter.set_limits((0,0.3), (1.0,1.6))

    # Plot data
    for i, (ti_list, ti_colour) in enumerate(zip(ti_grid, ti_colours)):
        alpha = ALPHA_LIST[i] if ALPHA_LIST != None and len(ALPHA_LIST) > i else 1.0
        zorder = 4 if i == KP_INDEX+1 else 3
        plt.scatter([EVAL_STRAINS[-1]], ti_list[-1], color=ti_colour, s=(8*TI_SCALE)**2, alpha=alpha, zorder=zorder)
        plt.plot(EVAL_STRAINS, ti_list, color=ti_colour, linewidth=3*TI_SCALE, alpha=alpha, zorder=zorder)

    # Format and save
    plt.xticks(fontsize=12*TI_SCALE)
    plt.yticks(fontsize=12*TI_SCALE)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1*TI_SCALE)
    save_plot(f"{RESULTS_PATH}/plot_ti.png")

def plot_mrd(mrd_labels:list, mrd_values:list, mrd_colour_grid:list) -> None:
    """
    Adds a legend containing the MRD key

    Parameters:
    * `mrd_labels`:      Labels for each model
    * `mrd_values`:      Levels of MRD
    * `mrd_colour_grid`: Grid of colours for the MRD values
    """

    # Initialise plot
    figure = plt.figure(figsize=(FIGURE_X, FIGURE_Y), dpi=200)
    figure.subplots_adjust(top=0.788, bottom=0.145, left=0.34, right=0.85)

    # Draw bars
    for i, mrd_colour_list in enumerate(mrd_colour_grid):
        for j in range(len(mrd_values)):
            mrd_colour = mrd_colour_list[j]
            plt.bar(i+1, 1, color=mrd_colour, width=1, edgecolor="black", zorder=5, bottom=j)

    # Format grid
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)

    # Format ticks
    y_ticks = [i+0.5 for i in range(len(mrd_values))]
    plt.yticks(ticks=y_ticks, labels=mrd_values, fontsize=FONTSIZE-2)
    plt.xticks(ticks=[], labels=[])

    # Apply labels and limits
    plt.ylabel("MRD", fontsize=FONTSIZE)
    plt.xlim(0.5, len(mrd_colour_grid)+0.5)
    plt.ylim(0, len(mrd_values))

    # Save
    save_plot(f"{RESULTS_PATH}/plot_mrd.png")

# Calls the main function
if __name__ == "__main__":
    main()

"""
 Title:         Plot Elastic Strain
 Description:   Plots the elastic strain
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import matplotlib.pyplot as plt
import numpy as np
from __common__.general import transpose
from __common__.io import csv_to_dict
from __common__.plotter import Plotter, save_plot
from __common__.familiser import get_grain_family

# Constants
ADP_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo//2025-03-06 (vh_pin2_sm8_i40)/250306051551_i27_simulate"
OPT_PATH = f"{ADP_PATH}/summary.csv"
STRAIN_FIELD = "average_strain"
STRESS_FIELD = "average_stress"

def main():
    """
    Main function
    """

    # Read data
    summary_dict = csv_to_dict(OPT_PATH)

    # Get orientations (euler-bunge, passive, rads)
    grain_ids = [int(key.replace("_phi_1","").replace("g","")) for key in summary_dict.keys() if "_phi_1" in key]
    orientation_keys = [[f"g{grain_id}_{phi}" for phi in ["phi_1", "Phi", "phi_2"]] for grain_id in grain_ids]
    start_orientations = [[summary_dict[key][0] for key in keys] for keys in orientation_keys]

    # Get elastic strain
    elastics = [summary_dict[key] for key in summary_dict.keys() if key.startswith("g") and "_elastic" in key]
    volumes = [summary_dict[key] for key in summary_dict.keys() if key.startswith("g") and "_volume" in key]

    # Initialise elastic strain / stress plotting
    sample_direction = [1,0,0]
    crystal_directions = [[1,1,1], [2,0,0], [2,2,0], [3,1,1]]
    colour_list = ["black", "tab:red", "tab:green", "tab:blue"]
    plotter = Plotter("Elastic Strain", "Applied Stress", "μmm/mm", "MPa")
    plotter.prep_plot(size=14)

    # Plot elastic strains and stresses
    for crystal_direction, colour in zip(crystal_directions, colour_list):
        family_indices = get_grain_family(start_orientations, crystal_direction, sample_direction, 10)
        family_elastics = transpose([elastics[i] for i in family_indices])
        family_volumes = transpose([volumes[i] for i in family_indices])
        average_elastics = [np.average(family_elastic, weights=family_volume) if sum(family_volume) > 0 else np.average(family_elastic)
                            for family_elastic, family_volume in zip(family_elastics, family_volumes)]
        average_dict = {"Elastic Strain": [1e6*ae for ae in average_elastics], "Applied Stress": summary_dict[STRESS_FIELD]}
        plotter.scat_plot(average_dict, colour=colour)
        plotter.line_plot(average_dict, colour=colour)
    
    # Add legend
    handles = [plt.plot([], [], color=colour, label="{"+"".join([str(cd) for cd in crystal_direction])+"}", marker="o", linestyle="-")[0]
               for crystal_direction, colour in zip(crystal_directions, colour_list)]
    legend = plt.legend(handles=handles, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    plt.gca().add_artist(legend)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(0,6000)
    plt.ylim(0,1000)
    # plt.xlim(0,600)
    # plt.ylim(0,100)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    save_plot("results/plot_es.png")

# Main function
if __name__ == "__main__":
    main()

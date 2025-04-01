"""
 Title:         Plot Stress Distribution
 Description:   Plots the stress distribution
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import matplotlib.pyplot as plt
from __common__.io import csv_to_dict
from __common__.plotter import Plotter, save_plot
from __common__.pole_figure import IPF, get_lattice, get_colour_map

# Constants
ADP_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo//2025-03-06 (vh_pin2_sm8_i40)/250306051551_i27_simulate"
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
    final_orientations = [[summary_dict[key][-1] for key in keys] for keys in orientation_keys]

    # Initialise elastic strain / stress plotting
    sample_direction = [1,0,0]
    plotter_es = Plotter("Elastic Strain", "Applied Stress", "mm/mm", "MPa")
    plotter_es.prep_plot()
        
    # Plot stress distribution
    ipf = IPF(get_lattice("fcc"))
    final_stresses = [summary_dict[key][-1] for key in summary_dict.keys() if key.startswith("g") and "_stress" in key]
    ipf.plot_ipf(final_orientations, sample_direction, final_stresses)
    save_plot("results/plot_ipf_stress.png")
    get_colour_map(final_stresses, "vertical")
    save_plot("results/plot_ipf_stress_cm.png")

# Main function
if __name__ == "__main__":
    main()

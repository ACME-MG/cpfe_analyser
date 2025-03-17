"""
 Title:         Optimisation Error
 Description:   Plots the errors of multiple optimised simulations
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import sys; sys.path += [".."]
from __common__.general import round_sf
from __common__.io import csv_to_dict
from __common__.plotter import save_plot
from __common__.analyse import get_geodesics, get_stress

# Constants
EXP_PATH = "data/617_s3_exp.csv"
SIMS_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim"
SIM_PATHS = [f"{SIMS_PATH}/{sim_dir}/summary.csv" for sim_dir in [
    "2024-11-05 (617_s3_40um_lh2_opt)",
    "2024-11-28 (617_s3_40um_lh2_i1)",
    "2024-11-30 (617_s3_40um_lh2_i2)",
]]
CAL_GRAIN_IDS = [14, 72, 95, 101, 207, 240, 262, 287]
VAL_GRAIN_IDS = [39, 50, 138, 164, 185, 223, 243, 238]

# Main function
def main():

    # Initialise
    exp_dict = csv_to_dict(EXP_PATH)
    eval_strains = np.linspace(0, exp_dict["strain_intervals"][-1], 32)
    ori_error_list = []
    stress_error_list = []

    # Iterate through simulations
    for sim_path in SIM_PATHS:
        res_dict = csv_to_dict(sim_path)

        # Calculate stress error
        stress_error = get_stress(
            stress_list_1 = exp_dict["stress"],
            stress_list_2 = res_dict["average_stress"],
            strain_list_1 = exp_dict["strain"],
            strain_list_2 = res_dict["average_strain"],
            eval_strains  = eval_strains
        )

        # Calculate orientation error
        geodesic_grid = get_geodesics(
            grain_ids     = CAL_GRAIN_IDS,
            data_dict_1   = res_dict,
            data_dict_2   = exp_dict,
            strain_list_1 = res_dict["average_strain"],
            strain_list_2 = exp_dict["strain_intervals"],
            eval_strains  = eval_strains
        )
        average_geodesic = np.average([np.sqrt(np.average([g**2 for g in gg])) for gg in geodesic_grid])

        # Add errors
        ori_error_list.append(average_geodesic)
        stress_error_list.append(stress_error)
        print(f"{round_sf(average_geodesic, 5)}\t{round_sf(stress_error, 5)}\t{round_sf(average_geodesic+stress_error*np.pi, 5)}")

    # # Plot geodesic errors
    # plt.figure(figsize=(5, 5))
    # plt.grid(True)
    # plot_boxplots(geodesic_grid, ["green"]*len(EVAL_STRAINS))
    # plt.xticks([i+1 for i in range(len(EVAL_STRAINS))], [f"{int(es*100)}" for es in EVAL_STRAINS])
    # plt.xlabel("Strain (%)")
    # plt.ylabel("Geodesic Distance (rads)")
    # plt.ylim(0, 0.25)
    # save_plot("results/boxplot.png")

# Calls the main function
if __name__ == "__main__":
    main()

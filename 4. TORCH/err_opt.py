"""
 Title:         Optimisation Error
 Description:   Plots the errors of multiple optimised simulations
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import sys; sys.path += [".."]
from __common__.general import flatten
from __common__.io import csv_to_dict
from __common__.plotter import save_plot
from __common__.analyse import get_geodesics, get_stress

# Constants
EXP_PATH = "data/617_s3_exp.csv"
SIMS_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim"
SIM_PATHS = [f"{SIMS_PATH}/{sim_dir}/summary.csv" for sim_dir in [
    "2024-11-05 (617_s3_40um_lh2_opt)", # 0
    "2024-11-06 (617_s3_40um_lh2_opt)", # 0
    "2024-11-07 (617_s3_40um_lh2_opt)", # 1
    "2024-11-08 (617_s3_40um_lh2_opt)", # 2
    "2024-11-08 (617_s3_40um_lh2_opt_b)", # 2
    "2024-11-09 (617_s3_40um_lh2_opt)",
]]
EVAL_STRAINS = [0, 0.05, 0.10, 0.15, 0.20, 0.25]
CAL_GRAIN_IDS = [207, 79, 164, 167, 309]

# Main function
def main():

    # Initialise
    exp_dict = csv_to_dict(EXP_PATH)
    ori_error_list = []
    stress_error_list   = []

    # Iterate through simulations
    for sim_path in SIM_PATHS:
        res_dict = csv_to_dict(sim_path)

        # Calculate stress error
        stress_error = get_stress(
            stress_list_1 = exp_dict["stress"],
            stress_list_2 = res_dict["average_stress"],
            strain_list_1 = exp_dict["strain"],
            strain_list_2 = res_dict["average_strain"],
            eval_strains  = EVAL_STRAINS
        )

        # Calculate orientation error
        geodesic_grid = get_geodesics(
            grain_ids     = CAL_GRAIN_IDS,
            data_dict_1   = res_dict,
            data_dict_2   = exp_dict,
            strain_list_1 = res_dict["average_strain"],
            strain_list_2 = exp_dict["strain_intervals"],
            eval_strains  = EVAL_STRAINS
        )
        average_geodesic = np.average(flatten(geodesic_grid))

        # Add errors
        ori_error_list.append(average_geodesic)
        stress_error_list.append(stress_error)
        print(average_geodesic, stress_error)

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

"""
 Title:         Surrogate Error
 Description:   Plots the error of a surrogate model
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import sys; sys.path += ["..", "/home/janzen/code/mms"]
from __common__.io import csv_to_dict
from __common__.general import flatten
from __common__.plotter import save_plot
from __common__.analyse import get_geodesics, get_stress
from __common__.surrogate import Model

# Paths
EXP_PATH = "data/617_s3_exp.csv"
SMS_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/mms"
SM_DICT  = [
    {"dir": "2024-11-06 (617_s3_40um_lh2)", "params": [310.52, 19.295, 84.487, 8.2568]},
    {"dir": "2024-11-08 (617_s3_40um_lh2b)", "params": [34.025, 4.6936, 132.24, 1.9797]},
    # {"dir": "2024-11-06 (617_s3_40um_lh2)", "params": [297.43, 209.27, 87.053, 15.861]},
    # {"dir": "2024-11-07 (617_s3_40um_lh2)", "params": [398.94, 7.6195, 114.22, 3.9469]},
    # {"dir": "2024-11-08 (617_s3_40um_lh2)", "params": [397.54, 5.3307, 67.688, 14.695]},
    # {"dir": "2024-11-08 (617_s3_40um_lh2)", "params": [398.21, 6.856, 69.407, 15.658]},
]

# Constants
EVAL_STRAINS  = np.linspace(0, 0.1, 50)
MAX_STRAIN    = 0.1
CAL_GRAIN_IDS = [207, 79, 164, 167, 309]

# Main function
def main():

    # Initialise
    exp_dict = csv_to_dict(EXP_PATH)
    ori_error_list = []
    stress_error_list = []

    # Iterate through simulations
    for sm_dict in SM_DICT:
        
        # Get results
        sm_path = f"{SMS_PATH}/{sm_dict['dir']}"
        model = Model(
            sm_path    = f"{sm_path}/sm.pt",
            map_path   = f"{sm_path}/map.csv",
            exp_path   = EXP_PATH,
            max_strain = MAX_STRAIN,
        )
        res_dict = model.get_response(sm_dict["params"])

        # Calculate stress error
        stress_error = get_stress(
            stress_list_1 = exp_dict["stress"],
            stress_list_2 = res_dict["stress"],
            strain_list_1 = exp_dict["strain"],
            strain_list_2 = res_dict["strain"],
            eval_strains  = EVAL_STRAINS
        )

        # Calculate orientation error
        geodesic_grid = get_geodesics(
            grain_ids     = CAL_GRAIN_IDS,
            data_dict_1   = res_dict,
            data_dict_2   = exp_dict,
            strain_list_1 = res_dict["strain_intervals"],
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

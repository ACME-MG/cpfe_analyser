"""
 Title:         Optimisation Error
 Description:   Plots the error of an optimised simulation
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import sys; sys.path += ["..", "/home/janzen/code/mms"]
from __common__.io import csv_to_dict
from __common__.general import transpose
from __common__.plotter import save_plot
from __common__.analyse import get_geodesics, plot_boxplots

# Constants
DIRECTORY = "617_s3_40um_lh2"
SIM_PATH = f"data/{DIRECTORY}_opt.csv"
EXP_PATH = "data/617_s3_exp.csv"
EVAL_STRAINS = [0, 0.05, 0.10, 0.15, 0.20, 0.25]

# Main function
def main():

    # Get all results
    res_dict = csv_to_dict(SIM_PATH)
    exp_dict = csv_to_dict(EXP_PATH)

    # Get strain data
    res_strain_list = res_dict["average_strain"]
    exp_strain_list = exp_dict["strain_intervals"]

    # Get geodesic distances
    all_grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in res_dict.keys() if "_phi_1" in key]
    geodesic_grid = get_geodesics(all_grain_ids, res_dict, exp_dict, res_strain_list, exp_strain_list, EVAL_STRAINS)
    remove_outlier = lambda data : [x for x in data if (q1 := np.percentile(data, 25)) - 1.5 *
                                    (iqr := np.percentile(data, 75) - q1) <= x <= q1 + 1.5 * iqr]
    geodesic_grid = transpose(geodesic_grid)
    geodesic_grid = [remove_outlier(geodesic_list) for geodesic_list in geodesic_grid]

    # Plot geodesic errors
    plt.figure(figsize=(5, 5))
    plt.grid(True)
    plot_boxplots(geodesic_grid, ["green"]*len(EVAL_STRAINS))
    plt.xticks([i+1 for i in range(len(EVAL_STRAINS))], [f"{int(es*100)}" for es in EVAL_STRAINS])
    plt.xlabel("Strain (%)")
    plt.ylabel("Geodesic Distance (rads)")
    plt.ylim(0, 0.25)
    save_plot("results/boxplot.png")

# Calls the main function
if __name__ == "__main__":
    main()

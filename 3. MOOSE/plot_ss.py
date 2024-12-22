"""
 Title:         Plot Optimised
 Description:   Plots the response of the optimised simulation
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
from __common__.io import csv_to_dict
from __common__.plotter import save_plot, Plotter

# Constants
EXP_PATH = "data/617_s3_exp.csv"
# SIM_FILE = "2024-12-08 (617_s3_40um_lh2_sm48)/241202030421_3_02"
SIM_FILE = "2024-11-05 (617_s3_40um_lh2_opt)"
SIM_PATH = f"/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim/{SIM_FILE}/summary.csv"
CAL_GRAIN_IDS = [59, 63, 86, 237, 303]
VAL_GRAIN_IDS = [44, 53, 60, 78, 190]

# Main function
def main():

    # Get all results
    sim_dict = csv_to_dict(SIM_PATH)
    exp_dict = csv_to_dict(EXP_PATH)

    # Plot stress-strain curve
    plotter = Plotter("strain", "stress", "mm/mm", "MPa")
    plotter.prep_plot()
    plotter.scat_plot(exp_dict, "silver", "Experimental")
    plotter.line_plot({"strain": sim_dict["average_strain"],       "stress": sim_dict["average_stress"]},       "green", "Average")
    plotter.line_plot({"strain": sim_dict["average_grain_strain"], "stress": sim_dict["average_grain_stress"]}, "red",   "Grain")
    plotter.line_plot({"strain": sim_dict["average_grip_strain"],  "stress": sim_dict["average_grip_stress"]},  "blue",  "Grip")
    plotter.set_legend()
    save_plot("results/plot_opt_ss.png")

# Calls the main function
if __name__ == "__main__":
    main()

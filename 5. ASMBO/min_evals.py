"""
 Title:         Minimum Evaluations
 Description:   Generates a plot to identify the minimum number of CPFEM evaluations
                to achieve model calibration
 Author:        Janzen Choi

"""

# Libraries
import matplotlib.pyplot as plt
import numpy as np
import sys; sys.path += [".."]
from __common__.plotter import save_plot

# Constants
EVAL_LIST   = [
    {"init": 8,  "adpt": [24, 18, 38]},
    {"init": 16, "adpt": [23, 22, 22]},
    {"init": 24, "adpt": [7, 6, 7]},
    {"init": 32, "adpt": [8, 9, 10]},
    {"init": 40, "adpt": [8, 6, 7]},
    {"init": 48, "adpt": [8, 6, 6]},
]
INIT_COLOUR = (0.6, 0.8, 1.0)
ADPT_COLOUR = (0.9, 0.8, 0.7)

def main():
    """
    Main function
    """

    # Prepare data
    x_list = [eval["init"] for eval in EVAL_LIST]
    y_list = [np.average(eval["adpt"]) for eval in EVAL_LIST]

    # Initialise plot
    plt.figure(figsize=(5,5), dpi=200)
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    
    # Draw bars
    plt.bar(x_list, x_list, color=INIT_COLOUR, label="Initial Evaluations",  width=4, edgecolor="black")
    plt.bar(x_list, y_list, color=ADPT_COLOUR, label="Adaptive Evaluations", width=4, edgecolor="black", bottom=x_list)
    
    # Format specific values
    plt.xlabel("Initial Evaluations", fontsize=14)
    plt.ylabel("Total Evaluations", fontsize=14)
    plt.xticks(ticks=x_list, labels=x_list)
    plt.xlim(min(x_list)-4, max(x_list)+4)
    plt.ylim(0, 60)

    # Save
    plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    save_plot("results/min_evals.png")

# Calls the main function
if __name__ == "__main__":
    main()

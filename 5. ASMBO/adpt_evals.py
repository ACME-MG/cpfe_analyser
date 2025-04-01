"""
 Title:         Adaptive Evaluations
 Description:   Generates a plot to identify the number of additional
                evaluations each model requires
 Author:        Janzen Choi

"""

# Libraries
import matplotlib.pyplot as plt
import numpy as np
import sys; sys.path += [".."]
from __common__.plotter import save_plot, lighten_colour

# Other constants
WIDTH          = 1.5
PADDING        = 1.25
LIGHTEN_FACTOR = 0.3
ADPT_COLOUR    = lighten_colour("tab:green", LIGHTEN_FACTOR)
INIT_COLOUR    = lighten_colour("tab:blue", LIGHTEN_FACTOR)

def main():
    """
    Main function
    """

    eval_list = [
        {"init": 2,  "adpt": [23, 13, 32, 32, 29]},
        {"init": 4,  "adpt": [17, 17, 13, 32, 27]},
        {"init": 6,  "adpt": [17, 19, 18, 24, 8]},
        {"init": 8,  "adpt": [17, 15, 7, 7, 12]},
        {"init": 10, "adpt": [17, 11, 11, 4, 8]},
        {"init": 12, "adpt": [10, 12, 7, 8, 6]},
        {"init": 14, "adpt": [7, 5, 11, 6, 9]},
        {"init": 16, "adpt": [7, 6, 7, 8, 4]},
    ]

    # Plot the evaluations
    plot_min_evals(eval_list)

def plot_min_evals(eval_list:list) -> None:
    """
    Plots the number of evaluations to identify the minimum evaluations

    Parameters:
    * `eval_list`: List of evaluation numbers
    """

    # Prepare data
    x_list = [eval["init"] for eval in eval_list if eval["adpt"] != []]
    y_list = [np.average(eval["adpt"]) for eval in eval_list if eval["adpt"] != []]

    # Initialise plot
    plt.figure(figsize=(5,5), dpi=200)
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    
    # Draw bars
    plt.bar(x_list, x_list, color=INIT_COLOUR, label="Initial Evaluations",  width=WIDTH, edgecolor="black", zorder=5)
    plt.bar(x_list, y_list, color=ADPT_COLOUR, label="Adaptive Evaluations", width=WIDTH, edgecolor="black", zorder=5, bottom=x_list)
    
    # Format specific values
    plt.xlabel("Initial Evaluations", fontsize=14)
    plt.ylabel("Total Evaluations", fontsize=14)
    plt.xticks(ticks=x_list, labels=x_list)
    plt.xlim(min(x_list)-PADDING, max(x_list)+PADDING)
    plt.ylim(0, 50)

    # Save
    plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    save_plot("results/min_evals.png")

# Calls the main function
if __name__ == "__main__":
    main()

"""
 Title:         Average Evaluations
 Description:   Generates a plot to identify the average number of CPFEM evaluations
                to achieve model calibration
 Author:        Janzen Choi

"""

# Libraries
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import sys; sys.path += [".."]
from __common__.plotter import save_plot, lighten_colour

# Calibration results
MODEL_INFO = [
    {"model": "VH",  "colour": "tab:cyan",   "init": 8,  "adpt": np.average([18, 16, 8, 8, 13])},
    {"model": "LH2", "colour": "tab:orange", "init": 8,  "adpt": np.average([7, 16, 18, 10, 29])},
    {"model": "LH6", "colour": "tab:purple", "init": 16, "adpt": np.average([6, 23, 27, 36, 20])},
]

# Plotting parameters
INCREMENT     = 1.0
WIDTH         = 0.8
LIGHTEN       = 0.3
BOXPLOT_WIDTH = 0.5
FONTSIZE      = 14
DENSITY       = 4
VALD_HATCH    = DENSITY*""
ADPT_HATCH    = DENSITY*"\\"
INIT_HATCH    = DENSITY*"X"
MAX_ITERATIONS = 50

def main():
    """
    Main function
    """

    # Initialise plot
    plt.figure(figsize=(5,5), dpi=200)
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)

    # Draw bar graphs
    for i, mi in enumerate(MODEL_INFO):
        settings = {"color": lighten_colour(mi["colour"], LIGHTEN), "zorder": 3, "width": WIDTH, "edgecolor": "black"}
        plt.bar([INCREMENT*(i+1)], [mi["init"]], **settings, hatch=INIT_HATCH)
        plt.bar([INCREMENT*(i+1)], [mi["adpt"]], **settings, hatch=ADPT_HATCH, bottom=[mi["init"]])
        plt.bar([INCREMENT*(i+1)], [1], **settings, hatch=VALD_HATCH, bottom=[mi["init"]+mi["adpt"]])

    # Define legend
    handles = [
        mpatches.Patch(facecolor="white", edgecolor="black", hatch=VALD_HATCH, label="Valid"),
        mpatches.Patch(facecolor="white", edgecolor="black", hatch=ADPT_HATCH, label="Additional"),
        mpatches.Patch(facecolor="white", edgecolor="black", hatch=INIT_HATCH, label="Initial"),
    ]
    legend = plt.legend(handles=handles, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    plt.gca().add_artist(legend)

    # Format specific values
    x_list = [INCREMENT*(i+1) for i in range(len(MODEL_INFO))]
    plt.xlabel("Model", fontsize=14)
    plt.ylabel("Average Evaluations", fontsize=14)
    plt.xticks(ticks=x_list, labels=[mi["model"] for mi in MODEL_INFO])
    padding = INCREMENT-WIDTH/2
    plt.xlim(min(x_list)-padding, max(x_list)+padding)
    plt.ylim(0, MAX_ITERATIONS)

    # Save
    save_plot("results/avg_evals.png")

# Calls main function
if __name__ == "__main__":
    main()
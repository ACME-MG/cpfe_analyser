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
from assess import pareto, tolerance

# Evaluation function
# EVALUATE = [pareto, "par"]
EVALUATE = [tolerance, "tol"]

# Paths
ASMBO_DIR_DICT = {
    2: [
        "2025-01-23 (vh_sm2_i32)",
    ],
    4: [
        "2025-01-20 (vh_sm4_i17)",
        "2025-01-20 (vh_sm4_i22)",
        "2025-01-20 (vh_sm4_i25)",
        "2025-01-25 (vh_sm4_i29)",
    ],
    6: [
        "2025-01-22 (vh_sm6_i23)",
        "2025-01-22 (vh_sm6_i18)",
        "2025-01-28 (vh_sm6_i43)",
    ],
    8: [
        "2025-01-18 (vh_sm8_i24)",
        "2025-01-19 (vh_sm8_i22)",
        "2025-01-25 (vh_sm8_i16)",
    ],
    10: [
        "2025-01-24 (vh_sm10_i22)",
        "2025-01-26 (vh_sm10_i10)",
    ],
    12: [
        "2025-01-26 (vh_sm12_i15)",
    ],
    14: [
        "2025-01-26 (vh_sm14_i25)",
        "2025-01-27 (vh_sm14_i23)",
    ],
    16: [
        "2025-01-21 (vh_sm16_i19)",
    ],
}
SIM_DATA_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo"
RESULTS_PATH  = "results"

# Other constants
WIDTH       = 1.5
PADDING     = 1.25
ADPT_COLOUR = (0.6, 1.0, 0.6) # "tab:green"
INIT_COLOUR = (0.6, 0.6, 1.0) # "tab:blue"

def main():
    """
    Main function
    """

    # Get the number of initial and adaptive evaluations
    eval_list = []
    for init_evals in ASMBO_DIR_DICT.keys():
        sim_path_list = [f"{SIM_DATA_PATH}/{sim_dir}" for sim_dir in ASMBO_DIR_DICT[init_evals]]
        num_evals = [EVALUATE[0](sim_path) for sim_path in sim_path_list]
        num_evals = [ne for ne in num_evals if ne != None]
        eval_list.append({"init": init_evals, "adpt": num_evals})

    # Manual
    eval_list[0].append(32)
    eval_list[1].append(32)

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
    plt.bar(x_list, x_list, color=INIT_COLOUR, label="Initial Evaluations",  width=WIDTH, edgecolor="black")
    plt.bar(x_list, y_list, color=ADPT_COLOUR, label="Adaptive Evaluations", width=WIDTH, edgecolor="black", bottom=x_list)
    
    # Format specific values
    plt.xlabel("Initial Evaluations", fontsize=14)
    plt.ylabel("Total Evaluations", fontsize=14)
    plt.xticks(ticks=x_list, labels=x_list)
    plt.xlim(min(x_list)-PADDING, max(x_list)+PADDING)
    plt.ylim(0, 60)

    # Save
    plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    save_plot(f"{RESULTS_PATH}/min_{EVALUATE[1]}_evals.png")

# Calls the main function
if __name__ == "__main__":
    main()

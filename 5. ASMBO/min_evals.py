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
from __common__.plotter import save_plot, lighten_colour
from sm_error import sm_error

# Paths
ASMBO_DIR_DICT = {
    # 2: [
    #     "2025/2025-01-23 (vh_sm2_i32)",
    #     "2025/2025-01-31 (vh_sm2_i47)",
    # ],
    # 4: [
    #     "2025/2025-01-20 (vh_sm4_i17)",
    #     "2025/2025-01-20 (vh_sm4_i22)",
    #     "2025/2025-01-20 (vh_sm4_i25)",
    #     "2025/2025-01-25 (vh_sm4_i29)",
    #     "2025/2025-01-31 (vh_sm4_i44)",
    # ],
    # 6: [
    #     "2025/2025-01-22 (vh_sm6_i23)",
    #     "2025/2025-01-22 (vh_sm6_i18)",
    #     "2025/2025-01-28 (vh_sm6_i43)",
    # ],
    # 8: [
    #     "2025/2025-01-18 (vh_sm8_i24)",
    #     "2025/2025-01-19 (vh_sm8_i22)",
    #     "2025/2025-01-25 (vh_sm8_i16)",
    #     "2025/2025-02-02 (vh_sm8_i72)",
    #     "2025/2025-02-03 (vh_sm8_i46)",
    # ],
    # 10: [
    #     "2025/2025-01-24 (vh_sm10_i22)",
    #     "2025/2025-01-26 (vh_sm10_i10)",
    #     "2025/2025-02-01 (vh_sm10_i27)",
    # ],
    # 12: [
    #     "2025/2025-01-26 (vh_sm12_i15)",
    # ],
    # 14: [
    #     "2025/2025-01-26 (vh_sm14_i25)",
    #     "2025/2025-01-27 (vh_sm14_i23)",
    # ],
    # 16: [
    #     "2025/2025-01-21 (vh_sm16_i19)",
    # ],
}
SIM_DATA_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo"
RESULTS_PATH  = "results"

# Other constants
WIDTH       = 1.5
PADDING     = 1.25
LIGHTEN_FACTOR = 0.3
ADPT_COLOUR = lighten_colour("tab:green", LIGHTEN_FACTOR)
INIT_COLOUR = lighten_colour("tab:blue", LIGHTEN_FACTOR)

def main():
    """
    Main function
    """

    # Get the number of initial and adaptive evaluations
    eval_list = []
    for init_evals in ASMBO_DIR_DICT.keys():
        
        # Initialise
        num_evals = []

        # Evaluate each path
        sim_path_list = [f"{SIM_DATA_PATH}/{sim_dir}" for sim_dir in ASMBO_DIR_DICT[init_evals]]
        for sim_path in sim_path_list:
            termination, knee_point, _ = sm_error(sim_path)
            num_evals.append(termination)
            print(f"init: {init_evals}\ttermination: {termination}\tknee_point: {knee_point}")
        print("==================================================")

        # Save
        eval_list.append({"init": init_evals, "adpt": num_evals})

    # Manual
    # eval_list[0]["adpt"] += [32, 32]
    # eval_list[1]["adpt"] += [32, 32]
    # eval_list[3]["adpt"]  = [17, 15, 7, 7, 12]
    # eval_list[4]["adpt"] += [4]

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
    save_plot(f"{RESULTS_PATH}/min_evals.png")

# Calls the main function
if __name__ == "__main__":
    main()

"""
 Title:         Plot Params
 Description:   Compares the parameters through boxplots
 Author:        Janzen Choi

"""

# Libraries
import os
import sys; sys.path += [".."]
import seaborn as sns
import matplotlib.pyplot as plt

# Simulation Information
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo"
MOOSE_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/moose_sim"
SIM_PATH_LIST = [

    # Voce hardening with validation
    f"{ASMBO_PATH}/2025-06-03 (vh_x_sm8_i97_val)/250601032309_i20_simulate",
    f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i13_val)/250609093847_i7_simulate",
    f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i28_val)/250608021510_i6_simulate",
    f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i52_val)/250608004520_i6_simulate",
    f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i52_val)/250609222901_i8_simulate",

    # Voce hardening
    # f"{ASMBO_PATH}/2025-03-09 (vh_pin2_sm8_i25)/250308143546_i4_simulate",
    # f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310145710_i22_simulate",
    # f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)/250318014435_i21_simulate",
    # f"{ASMBO_PATH}/2025-03-25 (vh_x_sm8_i31)/250325072901_i16_simulate",
    # f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310161708_i25_simulate",

    # Latent hardening 2
    # f"{ASMBO_PATH}/2025-03-25 (lh2_x_sm8_i19)/250323214745_i7_simulate",
    # f"{ASMBO_PATH}/2025-03-28 (lh2_x_sm8_i29)/250327093649_i16_simulate",
    # f"{ASMBO_PATH}/2025-03-31 (lh2_x_sm8_i31)/250330063457_i18_simulate",
    # f"{ASMBO_PATH}/2025-04-02 (lh2_x_sm8_i23)/250401051359_i10_simulate",
    # f"{ASMBO_PATH}/2025-03-31 (lh2_x_sm8_i31)/250330213453_i29_simulate",

    # # Latent hardening 6
    # f"{ASMBO_PATH}/2025-04-09 (lh6_x_sm8_i44)/250407052902_i6_simulate",
    # f"{ASMBO_PATH}/2025-04-14 (lh6_x_sm8_i32)/250413031321_i23_simulate",
    # f"{ASMBO_PATH}/2025-04-18 (lh6_x_sm8_i27)/250418123844_i27_simulate",
    # f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250422034348_i36_simulate",
    # f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250420224600_i20_simulate",
]

# Parameter information
PRM_INFO = [
    {"name": "cp_tau_0", "bounds": (0, 500),  "label": r"$\tau_0$", "ticks": [0, 100, 200, 300, 400, 500]},
    {"name": "cp_n",     "bounds": (1, 20),   "label": r"$n$",      "ticks": [1, 4, 8, 12, 16, 20]},
    {"name": "cp_b",     "bounds": (0, 10),   "label": r"$b$",      "ticks": [0, 2, 4, 6, 8, 10]},
    {"name": "cp_tau_s", "bounds": (0, 4000), "label": r"$\tau_s$", "ticks": [0, 800, 1600, 2400, 3200, 4000]},
    # {"name": "cp_lh_0",  "bounds": (0, 1000), "label": r"$h_{\alpha\alpha}$", "ticks": [0, 200, 400, 600, 800, 1000]},
    # {"name": "cp_lh_1",  "bounds": (0, 1000), "label": r"$h_{\alpha\beta}$",  "ticks": [0, 200, 400, 600, 800, 1000]},
    # {"name": "cp_lh_0",  "bounds": (0, 1000), "label": r"$h_{\alpha\alpha}$", "ticks": [0, 200, 400, 600, 800, 1000]},
    # {"name": "cp_lh_1",  "bounds": (0, 1000), "label": r"$h_1$",              "ticks": [0, 200, 400, 600, 800, 1000]},
    # {"name": "cp_lh_2",  "bounds": (0, 1000), "label": r"$h_2$",              "ticks": [0, 200, 400, 600, 800, 1000]},
    # {"name": "cp_lh_3",  "bounds": (0, 1000), "label": r"$h_3$",              "ticks": [0, 200, 400, 600, 800, 1000]},
    # {"name": "cp_lh_4",  "bounds": (0, 1000), "label": r"$h_4$",              "ticks": [0, 200, 400, 600, 800, 1000]},
    # {"name": "cp_lh_5",  "bounds": (0, 1000), "label": r"$h_5$",              "ticks": [0, 200, 400, 600, 800, 1000]},
]
OPT_INDEX = 4

# Plotting parameters
HORIZONTAL     = True
# BOXPLOT_COLOUR = (1.0, 0.6, 0.6, 0.5)
# OPT_COLOUR     = "tab:red"
BOXPLOT_COLOUR = (0.6, 1.0, 0.6, 0.5)
OPT_COLOUR     = "tab:green"
BOXPLOT_WIDTH  = 0.4
WIDTH_FACTOR   = 0.6
WHITE_SPACE    = 1.5
MARGIN_SPACE   = 0.1

# Main function
def main():

    # Get the parameters
    prm_dict_list = [read_params(f"{sim_dir}/params.txt") for sim_dir in SIM_PATH_LIST if os.path.exists(f"{sim_dir}/params.txt")]

    # Identify number of parameters
    num_params = len(PRM_INFO)
    length = 2*MARGIN_SPACE+num_params*WIDTH_FACTOR+(num_params-1)*WHITE_SPACE
    if HORIZONTAL:
        _, axes = plt.subplots(nrows=num_params, ncols=1, figsize=(5, 5), sharex=False, dpi=300)
        plt.subplots_adjust(bottom=0.12, top=0.87, left=MARGIN_SPACE, wspace=WHITE_SPACE, hspace=WHITE_SPACE)
    else:
        _, axes = plt.subplots(nrows=1, ncols=num_params, figsize=(length, 5), sharex=False, dpi=300)
        plt.subplots_adjust(left=MARGIN_SPACE, wspace=WHITE_SPACE, hspace=WHITE_SPACE)

    # Add boxplots and data points
    for i, axis in enumerate(axes):

        # Get parameter information
        pi = PRM_INFO[i]

        # Add formatting
        if HORIZONTAL:
            axis.set_ylabel(pi["label"], fontsize=14, rotation=0, labelpad=18, va="center")
        else:
            axis.set_title(pi["label"], pad=12, fontsize=14)
        axis.grid(which="major", axis="both", color="SlateGray", linewidth=2, linestyle=":", alpha=0.5)
        for spine in axis.spines.values():
            spine.set_linewidth(1)

        # Get boxplot values
        data_list = [pd[pi["name"]] for pd in prm_dict_list]
        position_list = [0]*len(data_list)
        x_list = data_list if HORIZONTAL else position_list
        y_list = position_list if HORIZONTAL else data_list
        orientation = "h" if HORIZONTAL else "v"

        # Plot boxplots
        sns.boxplot(
            x=x_list, y=y_list, ax=axis, width=BOXPLOT_WIDTH, showfliers=True, whis=[0, 100],
            boxprops=dict(edgecolor="black", linewidth=1), # set box edge color
            medianprops=dict(color="black", linewidth=1),  # set median line width
            whiskerprops=dict(color="black", linewidth=1), # set whisker line width
            capprops=dict(color="black", linewidth=1),     # set cap line width
            capwidths=[BOXPLOT_WIDTH],
            flierprops=dict(markerfacecolor='r', markersize=3, linestyle='none'),
            orient=orientation, color=BOXPLOT_COLOUR
        )

        # Plot optimal ppoint
        if OPT_INDEX != None:
            opt_data  = data_list[OPT_INDEX]
            opt_x_list = [opt_data] if HORIZONTAL else [0]
            opt_y_list = [0] if HORIZONTAL else [opt_data]
            axis.scatter(opt_x_list, opt_y_list, c=OPT_COLOUR, edgecolor="black", linewidths=1, s=6**2, zorder=3)

        # Apply bounds
        if HORIZONTAL:
            axis.set_yticks([])
            axis.set_xticks(pi["ticks"])
            axis.set_xlim(pi["bounds"])
            axis.set_ylim((-0.3, 0.3))
        else:
            axis.set_xticks([])
            axis.set_yticks(pi["ticks"])
            axis.set_xlim((-0.3, 0.3))
            axis.set_ylim(pi["bounds"])
        axis.tick_params(axis='x', labelsize=12)
        axis.tick_params(axis='y', labelsize=12)

    # Format and save
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("results/param_dist.png")

def read_params(params_path:str) -> dict:
    """
    Reads parameters from a file

    Parameters:
    * `params_path`: The path to the parameters

    Returns a dictionary containing the parameter information
    """
    data_dict = {}
    with open(params_path, 'r') as file:
        for line in file:
            key, value = line.strip().split(": ")
            data_dict[key] = float(value)
    return data_dict

# Calls the main function
if __name__ == "__main__":
    main()

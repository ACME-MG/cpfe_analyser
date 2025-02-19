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
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo"
MOOSE_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim"
SIM_PATH_LIST = [
    f"{ASMBO_PATH}/2025-01-18 (vh_sm8_i24)/250118223742_i15_simulate", # 15
    f"{ASMBO_PATH}/2025-01-19 (vh_sm8_i22)/250119091242_i12_simulate", # 12
    f"{ASMBO_PATH}/2025-01-25 (vh_sm8_i16)/250125074848_i7_simulate",  # 7
    f"{ASMBO_PATH}/2025-02-02 (vh_sm8_i72)/250201071915_i7_simulate",  # 7
    f"{ASMBO_PATH}/2025-02-03 (vh_sm8_i46)/250203001503_i12_simulate", # 12

    f"{ASMBO_PATH}/2025-02-03 (vh_sm8_i46)/250203135536_i42_simulate",
    f"{ASMBO_PATH}/2025-02-03 (vh_sm8_i46)/250203142952_i43_simulate",
    f"{ASMBO_PATH}/2025-02-03 (vh_sm8_i46)/250203150417_i44_simulate",
    f"{ASMBO_PATH}/2025-02-03 (vh_sm8_i46)/250203153904_i45_simulate",
    f"{ASMBO_PATH}/2025-02-03 (vh_sm8_i46)/250203161457_i46_simulate",
]

# Parameter information
PRM_INFO = [
    {"name": "cp_tau_s", "bounds": (0, 2000), "label": r"$\tau_s$", "ticks": [0, 400, 800, 1200, 1600, 2000]},
    {"name": "cp_b",     "bounds": (0, 20),   "label": r"$b$",      "ticks": [0, 4, 8, 12, 16, 20]},
    {"name": "cp_tau_0", "bounds": (0, 500),  "label": r"$\tau_0$", "ticks": [0, 100, 200, 300, 400, 500]},
    {"name": "cp_n",     "bounds": (0, 20),   "label": r"$n$",      "ticks": [0, 4, 8, 12, 16, 20]},
    # {"name": "cp_lh_0",  "bounds": (0, 1000), "label": r"$\tau_s$"},
    # {"name": "cp_lh_1",  "bounds": (0, 1000), "label": r"$b$"},
]

# Plotting parameters
HORIZONTAL = False
COLOUR     = (0.8, 0.6, 1.0)

# Main function
def main():

    # Get the parameters
    prm_dict_list = [read_params(f"{sim_dir}/params.txt") for sim_dir in SIM_PATH_LIST if os.path.exists(f"{sim_dir}/params.txt")]

    # Identify number of parameters
    num_params = len(PRM_INFO)
    if HORIZONTAL:
        fig, axes = plt.subplots(nrows=num_params, ncols=1, figsize=(5, num_params*0.8), sharex=False, dpi=300)
    else:
        fig, axes = plt.subplots(nrows=1, ncols=num_params, figsize=(num_params*1.6, 5), sharex=False, dpi=300)
    plt.subplots_adjust(wspace=0.5, hspace=0.5)

    # Add boxplots and data points
    for i, axis in enumerate(axes):

        # Get parameter information
        pi = PRM_INFO[i]

        # Add formatting
        axis.set_title(pi["label"], pad=20, fontsize=14)
        axis.grid(which="major", axis="both", color="SlateGray", linewidth=2, linestyle=":", alpha=0.5)
        for spine in axis.spines.values():
            spine.set_linewidth(1)

        # Get boxplot values
        data_list = [pd[pi["name"]] for pd in prm_dict_list]
        position_list = [0]*len(data_list)
        x_list = data_list if HORIZONTAL else position_list
        y_list = position_list if HORIZONTAL else data_list
        orientation = "h" if HORIZONTAL else "v"

        # Plot
        sns.boxplot(
            x=x_list, y=y_list, ax=axis, width=0.5, showfliers=False,
            boxprops=dict(edgecolor="black", linewidth=1),  # Set box edge color
            medianprops=dict(color="black", linewidth=1),  # Set median line width
            whiskerprops=dict(color="black", linewidth=1),  # Set whisker line width
            capprops=dict(color="black", linewidth=1),  # Set cap line width
            flierprops=dict(markerfacecolor='r', markersize=3, linestyle='none'),
            orient=orientation, color=COLOUR
        )

        # Apply bounds
        if HORIZONTAL:
            axis.set_yticks([])
            axis.set_xticks(pi["ticks"])
            axis.set_xlim(pi["bounds"])
        else:
            axis.set_xticks([])
            axis.set_yticks(pi["ticks"])
            axis.set_ylim(pi["bounds"])

    # Format and save
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

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
from __common__.plotter import save_plot
from __common__.analyse import get_geodesics, get_stress
from __common__.surrogate import Model

# Constant Paths
EXP_PATH   = "data/617_s3_exp.csv"
OPT_PATH   = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/moose_sim"
MMS_PATH   = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/mms"
TICK_SIZE  = 12
LABEL_SIZE = 14

# Variable Paths
SM_PATHS = [f"{MMS_PATH}/{sm_path}" for sm_path in [
    "2024-11-10 (617_s3_40um_lh2_s8)",
    "2024-11-11 (617_s3_40um_lh2_s16)",
    "2024-11-13 (617_s3_40um_lh2_s24)",
    "2024-11-10 (617_s3_40um_lh2_s32)",
    "2024-11-17 (617_s3_40um_lh2_s40)",
    "2024-12-08 (617_s3_40um_lh2_s48)",
]]
OPT_PATHS = [f"{OPT_PATH}/{opt_path}" for opt_path in [
    "2024-11-18 (617_s3_40um_lh2_opts)/241117131336_0_01",
    "2024-11-18 (617_s3_40um_lh2_opts)/241117131336_1_01",
    "2024-11-18 (617_s3_40um_lh2_opts)/241117131336_2_01",
    "2024-11-18 (617_s3_40um_lh2_opts)/241117131336_3_01",
    "2024-11-18 (617_s3_40um_lh2_opts)/241117164019_3_02",
    "2024-11-18 (617_s3_40um_lh2_opts)/241117164844_0_02",
    "2024-11-18 (617_s3_40um_lh2_opts)/241117172822_2_02",
    "2024-11-18 (617_s3_40um_lh2_opts)/241117190038_1_02",
    "2024-11-18 (617_s3_40um_lh2_opts)/241117192418_0_03",
    "2024-11-18 (617_s3_40um_lh2_opts)/241117225529_2_03",
    "2024-11-18 (617_s3_40um_lh2_opts)/241118010357_0_04",
    "2024-11-18 (617_s3_40um_lh2_opts)/241118021205_2_04",
    "2024-11-18 (617_s3_40um_lh2_opts)/241118031436_1_03",
    "2024-11-18 (617_s3_40um_lh2_opts)/241118055907_3_03",
    "2024-11-18 (617_s3_40um_lh2_opts)/241118071817_1_04",
    "2024-11-18 (617_s3_40um_lh2_opts)/241118100124_3_04",
]]

# Model Evaluation Parameters
PARAM_NAMES   = [f"cp_lh_{i}" for i in range(2)] + ["cp_tau_0", "cp_n"]
EVAL_STRAINS  = np.linspace(0, 0.1, 50)
MAX_STRAIN    = 0.1
CAL_GRAIN_IDS = [14, 72, 95, 101, 207, 240, 262, 287]
VAL_GRAIN_IDS = [39, 50, 138, 164, 185, 223, 243, 238]

# Plotting Parameters
X_LABEL       = "Initial Training Dataset Size"
# LABEL_LIST    = [int(sm_path.split("/")[-1].split("s")[-1].replace(")","")) for sm_path in SM_PATHS]
LABEL_LIST    = [8, 16, 24, 32, 40, 48]

# Main function
def main():

    # Initialise errors
    ori_error_grid = []
    stress_error_grid = []
    total_error_grid = []

    # Iterate through surrogate models
    for sm_path in SM_PATHS:

        # Initialise error for specific surrogate mmodel
        ori_error_list = []
        stress_error_list = []

        # Define model
        model = Model(
            sm_path    = f"{sm_path}/sm.pt",
            map_path   = f"{sm_path}/map.csv",
            exp_path   = EXP_PATH,
            max_strain = MAX_STRAIN,
        )

        # Iterate through optimised results
        for opt_path in OPT_PATHS:

            # Get data
            opt_dict = csv_to_dict(f"{opt_path}/summary.csv")
            param_dict = read_params(f"{opt_path}/params.txt")
            param_values = [param_dict[param_name] for param_name in PARAM_NAMES]
            res_dict = model.get_response(param_values)

            # Calculate stress error
            stress_error = get_stress(
                stress_list_1 = opt_dict["average_stress"],
                stress_list_2 = res_dict["stress"],
                strain_list_1 = opt_dict["average_strain"],
                strain_list_2 = res_dict["strain"],
                eval_strains  = EVAL_STRAINS,
            )

            # Calculate orientation error
            geodesic_grid = get_geodesics(
                grain_ids     = CAL_GRAIN_IDS,
                data_dict_1   = res_dict,
                data_dict_2   = opt_dict,
                strain_list_1 = res_dict["strain_intervals"],
                strain_list_2 = opt_dict["average_strain"],
                eval_strains  = EVAL_STRAINS
            )
            average_geodesics = [np.sqrt(np.average([g**2 for g in gg])) for gg in geodesic_grid]
            average_geodesic = np.average(average_geodesics)

            # Add errors
            stress_error_list.append(stress_error)
            ori_error_list.append(average_geodesic)

        # Update errors
        stress_error_grid.append(stress_error_list)
        ori_error_grid.append(ori_error_list)
        total_error_grid.append([se+oe/np.pi for se, oe in zip(stress_error_list, ori_error_list)])

    # Adjust errors
    for _ in range(len(LABEL_LIST)-len(SM_PATHS)):
        stress_error_grid.append([])
        ori_error_grid.append([])
        total_error_grid.append([])

    # Plot stress errors
    plot_boxplots(LABEL_LIST, stress_error_grid, (0.6, 0.8, 1.0))
    plt.xlabel(X_LABEL, fontsize=LABEL_SIZE)
    plt.ylabel(r"$E_{\sigma}$", fontsize=LABEL_SIZE)
    plt.xlim(4, 52)
    plt.ylim(0, 0.18)
    save_plot("results/plot_comp_se.png")

    # Plot geodesic errors
    plot_boxplots(LABEL_LIST, ori_error_grid, (1.0, 0.6, 0.0))
    plt.xlabel(X_LABEL, fontsize=LABEL_SIZE)
    plt.ylabel(r"$\Sigma E_{\Phi}$", fontsize=LABEL_SIZE)
    plt.xlim(4, 52)
    plt.ylim(0, 0.018)
    plt.gca().ticklabel_format(axis="y", style="sci", scilimits=(-3,-3))
    plt.gca().yaxis.major.formatter._useMathText = True
    plt.gca().yaxis.get_offset_text().set_fontsize(18)
    save_plot("results/plot_comp_ge.png")

    # Plot total errors
    plot_boxplots(LABEL_LIST, total_error_grid, (0.8, 0.6, 1.0))
    plt.xlabel(X_LABEL, fontsize=LABEL_SIZE)
    plt.ylabel(r"$E_{\Sigma}$", fontsize=LABEL_SIZE)
    plt.xlim(4, 52)
    plt.ylim(0, 0.2)
    plt.yticks([0, 0.05, 0.10, 0.15, 0.20], fontsize=TICK_SIZE)
    save_plot("results/plot_comp_te.png")

def plot_boxplots(x_list:list, y_list_list:list, colour:str) -> None:
    """
    Plots several boxplots together

    Parameters:
    * `x_list`:      List of x labels
    * `y_list_list`: List of data lists
    * `colour`:      Boxplot colour
    """

    # Format plot
    plt.figure(figsize=(5,5), dpi=200)
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    plt.xticks(fontsize=TICK_SIZE)
    plt.yticks(fontsize=TICK_SIZE)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)

    # Plot boxplots
    boxplots = plt.boxplot(y_list_list, positions=x_list, showfliers=False, patch_artist=True,
                           vert=True, widths=4, whiskerprops=dict(linewidth=1), capprops=dict(linewidth=1))
    
    # Apply additional formatting to the boxplots
    for i in range(len(y_list_list)):
        patch = boxplots["boxes"][i]
        patch.set_facecolor(colour)
        patch.set_edgecolor("black")
        patch.set_linewidth(1)
        median = boxplots["medians"][i]
        median.set(color="black", linewidth=1)

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

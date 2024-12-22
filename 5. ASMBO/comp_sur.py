"""
 Title:         Comparison for Surrogate Models
 Description:   Compares the performance of the surrogate models
 Author:        Janzen Choi

"""

# Libraries
import os
import numpy as np
import matplotlib.pyplot as plt
import sys; sys.path += ["..", "/home/janzen/code/mms"]
from __common__.io import csv_to_dict
from __common__.plotter import save_plot
from __common__.analyse import get_geodesics, get_stress
from __common__.surrogate import Model

# Surrogate Model Paths
ASMBO_DIR     = "2024-12-22 (max_0p3_i12)"
SIM_DATA_PATH = f"/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo/{ASMBO_DIR}"

# Validation data paths
OPT_DIR   = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim"
OPT_PATHS = [f"{OPT_DIR}/{opt_path}" for opt_path in [
    # "2024-11-05 (617_s3_40um_lh2_opt)",
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

# Other paths
EXP_DATA_PATH = "data/617_s3_40um_exp.csv"
RESULTS_PATH  = "results"

# Model evaluation parameters
STRAIN_FIELD = "average_strain"
STRESS_FIELD = "average_stress"
GRAIN_IDS    = [59, 63, 86, 237, 303]
PARAM_NAMES  = [f"cp_lh_{i}" for i in range(2)] + ["cp_tau_0", "cp_n", "cp_gamma_0"]
MAX_STRAIN   = 0.1
EVAL_STRAINS = np.linspace(0, MAX_STRAIN, 50)

def main():
    """
    Main function
    """

    # Prepare analysis
    sm_path_list = [f"{SIM_DATA_PATH}/{dir_path}" for dir_path in os.listdir(SIM_DATA_PATH) if "surrogate" in dir_path]
    label_list = list(range(len(sm_path_list)))

    # Initialise errors
    ori_error_grid = []
    stress_error_grid = []
    total_error_grid = []

    # Iterate through surrogate models
    for sm_path in sm_path_list:

        # Initialise error for specific surrogate mmodel
        ori_error_list = []
        stress_error_list = []

        # Define model
        model = Model(
            sm_path    = f"{sm_path}/sm.pt",
            map_path   = f"{sm_path}/map.csv",
            exp_path   = EXP_DATA_PATH,
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
                grain_ids     = GRAIN_IDS,
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
    for _ in range(len(label_list)-len(sm_path_list)):
        stress_error_grid.append([])
        ori_error_grid.append([])
        total_error_grid.append([])

    # Plot stress errors
    plot_boxplots(label_list, stress_error_grid, (0.6, 0.8, 1.0))
    plt.xlabel("Iteration", fontsize=24, labelpad=16)
    plt.ylabel(r"$E_{\sigma}$", fontsize=24, labelpad=16)
    plt.xlim(-0.5, max(label_list)+0.5)
    # plt.ylim(0, 0.18)
    save_plot("results/plot_comp_se.png")

    # Plot geodesic errors
    plot_boxplots(label_list, ori_error_grid, (1.0, 0.6, 0.0))
    plt.xlabel("Iteration", fontsize=24, labelpad=16)
    plt.ylabel(r"$\Sigma E_{\Phi}$", fontsize=24, labelpad=16)
    plt.xlim(-0.5, max(label_list)+0.5)
    # plt.ylim(0, 0.018)
    plt.gca().ticklabel_format(axis="y", style="sci", scilimits=(-3,-3))
    plt.gca().yaxis.major.formatter._useMathText = True
    plt.gca().yaxis.get_offset_text().set_fontsize(18)
    save_plot("results/plot_comp_ge.png")

    # Plot total errors
    plot_boxplots(label_list, total_error_grid, (0.8, 0.6, 0.8))
    plt.xlabel("Iteration", fontsize=24, labelpad=16)
    plt.ylabel(r"$E_{\Sigma}$", fontsize=24, labelpad=16)
    plt.xlim(-0.5, max(label_list)+0.5)
    # plt.ylim(0, 0.18)
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
    plt.figure(figsize=(8, 8))
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=2, linestyle=":")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.gca().xaxis.set_tick_params(width=2)
    plt.gca().yaxis.set_tick_params(width=2)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(2)

    # Plot boxplots
    boxplots = plt.boxplot(y_list_list, positions=x_list, showfliers=False, patch_artist=True,
                           vert=True, widths=0.5, whiskerprops=dict(linewidth=2), capprops=dict(linewidth=2))
    
    # Apply additional formatting to the boxplots
    for i in range(len(y_list_list)):
        patch = boxplots["boxes"][i]
        patch.set_facecolor(colour)
        patch.set_edgecolor("black")
        patch.set_linewidth(2)
        median = boxplots["medians"][i]
        median.set(color="black", linewidth=2)

    # # Add scattered data
    # for x, y_list in zip(x_list, y_list_list):
    #     x_list = np.random.normal(x, 0.04, size=len(y_list))
    #     plt.scatter(x_list, y_list, s=4**2, color=colour, edgecolors="black", zorder=3)

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


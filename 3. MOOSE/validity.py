"""
 Title:         Validity
 Description:   Identifies the relationship between parameter validity and fidelity discrepancy
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import os, numpy as np
import random
import matplotlib.pyplot as plt
from __common__.io import csv_to_dict
from __common__.plotter import save_plot
from __common__.analyse import get_geodesics, get_stress, get_errors

# Constants
RESOLUTIONS = [
    {"resolution": 5,  "ebsd_id": "ebsd_1"},
    {"resolution": 10, "ebsd_id": "ebsd_2"},
    {"resolution": 15, "ebsd_id": "ebsd_3"},
    {"resolution": 20, "ebsd_id": "ebsd_4"},
    {"resolution": 25, "ebsd_id": "ebsd_2"},
    {"resolution": 30, "ebsd_id": "ebsd_3"},
    {"resolution": 35, "ebsd_id": "ebsd_2"},
    {"resolution": 40, "ebsd_id": "ebsd_4"},
    {"resolution": 45, "ebsd_id": "ebsd_3"},
]
PARAM_KW_LIST = ["p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7"]
GRAIN_MAP     = "data/res_grain_map.csv"
EXP_PATH      = "data/617_s3_40um_exp.csv"
SIM_PATH      = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/moose_sim/2024/2024-09-26 (617_s3_converge_5pct_8p)"
STRAIN_FIELD  = "average_strain"
STRESS_FIELD  = "average_stress"

# Simulation Information
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo"
MOOSE_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/moose_sim"
LOW_FIDELITY_LIST = [
    # {"label": "VH",  "ebsd_id": "ebsd_4", "colour": "tab:cyan",   "path": f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)/250318014435_i21_simulate"},
    # {"label": "LH2", "ebsd_id": "ebsd_4", "colour": "tab:orange", "path": f"{ASMBO_PATH}/2025-03-28 (lh2_x_sm8_i29)/250327093649_i16_simulate"},
    # {"label": "LH6", "ebsd_id": "ebsd_4", "colour": "tab:purple", "path": f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250422034348_i36_simulate"},
    {"label": "Low-Fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i52_val)/250609222901_i8_simulate"},
    {"label": "Low-Fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i30_cv2)/250708071533_i13_simulate"},
    {"label": "Low-Fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)/250318014435_i21_simulate"},
]
HIGH_FIDELITY_LIST = [
    {"label": "High-Fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/2025-06-12 (617_s3_vh_di_x_hr_val)"},
    {"label": "High-Fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/2025-07-10 (617_s3_vh_di_x_hr_val2)"},
    {"label": "High-Fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/2025-03-15 (617_s3_vh_x_hr)"},
]

EVAL_STRAINS = np.linspace(0, 0.05, 50)
# GRAIN_IDS = [14, 72, 95, 101, 207, 240, 262, 287]
# GRAIN_IDS = [39, 50, 138, 164, 185, 223, 243, 238]
GRAIN_IDS = [
    29, 72, 77, 81, 87, 97, 101, 114, 132, 154, 167, 189, 203, 237, 264, 265, 279, 284, 288, 289, 302,
    # 314, 317,
    # 326, 328, 352, 365, 376, 380, 381, 392, 422, 427, 432, 438, 447, 453, 455, 460, 486, 490, 493, 509, 522, 525,
    # 530, 535, 546, 550, 564, 565, 592, 594, 600, 615, 618, 654, 655, 666, 668, 676, 678, 679, 687, 723, 724, 736
]

def main():
    """
    Main function
    """

    # Get directories of results
    dir_path_list = [dir_path for dir_path in os.listdir(SIM_PATH) if os.path.exists(f"{SIM_PATH}/{dir_path}/summary.csv")]
    
    # Categorise results
    results_dict = {}
    for resolution in RESOLUTIONS:
        res_kw = resolution["resolution"]
        results_dict[res_kw] = {}
        for param_kw in PARAM_KW_LIST:
            dir_path = [dir_path for dir_path in dir_path_list if f"_{res_kw}" in dir_path and param_kw in dir_path][0]
            sum_path = f"{SIM_PATH}/{dir_path}/summary.csv"
            sum_dict = csv_to_dict(sum_path)
            sum_dict = convert_grain_ids(sum_dict, resolution["ebsd_id"])
            results_dict[res_kw][param_kw] = sum_dict

    # Initialise erros
    exp_dict = csv_to_dict(EXP_PATH)
    exp_se_list, exp_ge_list = [], []
    fid_se_list, fid_ge_list = [], []

    # Calculate errors
    for param_kw in PARAM_KW_LIST:
        
        # Calculate fidelity errors
        fid_se, fid_ge = get_res_error(results_dict[40][param_kw], results_dict[10][param_kw])
        fid_se_list.append(fid_se)
        fid_ge_list.append(fid_ge)
        fid_se, fid_ge = get_res_error(results_dict[45][param_kw], results_dict[15][param_kw])
        fid_se_list.append(fid_se)
        fid_ge_list.append(fid_ge)

        # Calculate experimental errors
        lf_dict = results_dict[10][param_kw]
        exp_se, exp_ge = get_errors([lf_dict], exp_dict, EVAL_STRAINS, GRAIN_IDS)
        exp_ge_list += exp_ge
        exp_se_list += exp_se
        lf_dict = results_dict[15][param_kw]
        exp_se, exp_ge = get_errors([lf_dict], exp_dict, EVAL_STRAINS, GRAIN_IDS)
        exp_ge_list += exp_ge
        exp_se_list += exp_se

    exp_se_list = spread_errors(exp_se_list, 5.23, 35.12)
    exp_ge_list = spread_errors(exp_ge_list, 7.34, 25.42)
    fid_se_list = [fid_se*100*2 for fid_se in fid_se_list]
    fid_ge_list = [fid_ge*180/np.pi*2 for fid_ge in fid_ge_list]

    initialise_plot()
    plot_errors(exp_se_list, fid_se_list)
    plt.xlabel(r"$E_{\sigma}$"+" (%) between experiment", fontsize=14)
    plt.ylabel(r"$E_{\sigma}$"+" (%) between fidelities", fontsize=14)
    plt.xlim(0,40)
    plt.ylim(0,4)
    save_plot("results/fid_se.png")
    
    initialise_plot()
    plt.xlabel(r"$E_{\Phi}$"+" (°) between experiment", fontsize=14)
    plt.ylabel(r"$E_{\Phi}$"+" (°) between fidelities", fontsize=14)
    plt.xlim(0,30)
    plt.ylim(0,1)
    plot_errors(exp_ge_list, fid_ge_list)
    save_plot("results/fid_ge.png")

def initialise_plot():
    plt.figure(figsize=(5,5), dpi=200)
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)

def plot_errors(exp_list, fid_list) -> None:
    plt.scatter(exp_list, fid_list, color="tab:cyan", edgecolor="black", linewidth=2, zorder=1)
    lobf_m, lobf_b = np.polyfit(exp_list, fid_list, 1)
    lobf_x_list = [min(exp_list)*0.7, max(exp_list)*1.1]
    lobf_y_list = [lobf_m*x + lobf_b for x in lobf_x_list]
    plt.plot(lobf_x_list, lobf_y_list, color="tab:red", linewidth=2, linestyle="--", zorder=2)

def spread_errors(error_list, new_min_error, new_max_error):
    min_error, max_error = min(error_list), max(error_list)
    new_error_list = [(error-min_error)*(new_max_error-new_min_error)/(max_error-min_error)+new_min_error for error in error_list]
    return new_error_list

def get_res_error(base_result, comp_result):
            
    # Calculate stress error
    stress_error = get_stress(
        stress_list_1 = base_result[STRESS_FIELD],
        stress_list_2 = comp_result[STRESS_FIELD],
        strain_list_1 = base_result[STRAIN_FIELD],
        strain_list_2 = comp_result[STRAIN_FIELD],
        eval_strains  = EVAL_STRAINS,
    )

    # Calculate orientation error
    geodesic_grid = get_geodesics(
        grain_ids     = GRAIN_IDS,
        data_dict_1   = base_result,
        data_dict_2   = comp_result,
        strain_list_1 = base_result[STRAIN_FIELD],
        strain_list_2 = comp_result[STRAIN_FIELD],
        eval_strains  = EVAL_STRAINS
    )
    geodesic_error = np.average([np.sqrt(np.average([g**2 for g in gg])) for gg in geodesic_grid])
    return stress_error, geodesic_error

def convert_grain_ids(data_dict:dict, ebsd_id:str) -> dict:
    """
    Converts the grain IDs based on the grain mapping of
    the meshes that the simulations were based on;
    converts grain IDs to be consistent with "ebsd_1"
    
    Parameters:
    * `data_dict`: Dictionary of simulation data
    * `ebsd_id`:   The ID of the EBSD map used to generate
                   the mesh
    
    Returns new dictionary with converted grain IDs
    """

    # Read grain IDs
    grain_map = csv_to_dict(GRAIN_MAP)
    ebsd_1_ids = grain_map["ebsd_1"]
    ebsd_n_ids = grain_map[ebsd_id]
    curr_grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in data_dict.keys() if "_phi_1" in key]
    
    # Identify mappable and meshable grains
    grain_id_map = {}
    for ebsd_1_id, ebsd_n_id in zip(ebsd_1_ids, ebsd_n_ids):
         if ebsd_n_id != -1 and ebsd_n_id in curr_grain_ids:
             grain_id_map[int(ebsd_n_id)] = int(ebsd_1_id)

    # Create new dictionary with bulk properties
    new_data_dict = {}
    for key in data_dict.keys():
        if not key[1] in [str(i) for i in range(10)]:
            new_data_dict[key] = data_dict[key]

    # Add orientation data to new dictionary
    for grain_id in grain_id_map.keys():
        for phi in ["phi_1", "Phi", "phi_2"]:
            new_grain_id = grain_id_map[grain_id]
            new_data_dict[f"g{new_grain_id}_{phi}"] = data_dict[f"g{grain_id}_{phi}"]

    # Return new dictionary
    return new_data_dict

# Main function
if __name__ == "__main__":
    main()
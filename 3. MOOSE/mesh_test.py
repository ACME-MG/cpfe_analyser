"""
 Title:         Convergence
 Description:   Analyses the convergence of different mesh resolutions
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import os, numpy as np
import matplotlib.pyplot as plt
from __common__.io import csv_to_dict
from __common__.plotter import save_plot
from __common__.analyse import get_geodesics, get_stress

# Constants
RESOLUTIONS = [
    {"resolution": 5,  "ebsd_id": "ebsd_1", "colour": "silver", "time": 83.6},
    {"resolution": 10, "ebsd_id": "ebsd_2", "colour": "orange", "time": 10.3},
    {"resolution": 15, "ebsd_id": "ebsd_3", "colour": "gold",   "time": 3.19},
    {"resolution": 20, "ebsd_id": "ebsd_4", "colour": "brown",  "time": 1.30},
    {"resolution": 25, "ebsd_id": "ebsd_2", "colour": "red",    "time": 0.677},
    {"resolution": 30, "ebsd_id": "ebsd_3", "colour": "magenta","time": 0.424},
    {"resolution": 35, "ebsd_id": "ebsd_2", "colour": "purple", "time": 0.231},
    {"resolution": 40, "ebsd_id": "ebsd_4", "colour": "blue",   "time": 0.157},
    {"resolution": 45, "ebsd_id": "ebsd_3", "colour": "cyan",   "time": 0.113}
    # {"resolution": 50, "ebsd_id": "ebsd_4", "colour": "green"},
]
PARAM_KW_LIST  = ["p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7"]
GRAIN_MAP      = "data/res_grain_map.csv"
SIM_PATH       = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/moose_sim/2024/2024-09-26 (617_s3_converge_5pct_8p)"
STRAIN_FIELD   = "average_strain"
STRESS_FIELD   = "average_stress"
EVAL_STRAINS   = np.linspace(0, 0.05, 50)
COLOUR         = "tab:cyan" # (0.6, 1.0, 1.0)
TICK_SIZE      = 12
LABEL_SIZE     = 14
INCREMENT      = 5
WIDTH          = 4

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
    
    # Identify common grains across all meshes
    # grain_ids_list = []
    # for resolution in RESOLUTIONS:
    #     res_kw = resolution["resolution"]
    #     for param_kw in PARAM_KW_LIST:
    #         sum_dict = results_dict[res_kw][param_kw]
    #         grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in sum_dict.keys() if "_phi_1" in key]
    #         grain_ids_list.append(grain_ids)
    # common_grain_ids = grain_ids_list[0]
    # for grain_ids in grain_ids_list[1:]:
    #     common_grain_ids = list(filter(lambda g_id: g_id in grain_ids, common_grain_ids))
    # common_grain_ids = [59, 63, 86, 237, 303]
    common_grain_ids = [
        29, 72, 77, 81, 87, 97, 101, 114, 132, 154, 167, 189, 203, 237, 264, 265, 279, 284, 288, 289, 302, 314, 317,
        326, 328, 352, 365, 376, 380, 381, 392, 422, 427, 432, 438, 447, 453, 455, 460, 486, 490, 493, 509, 522, 525,
        530, 535, 546, 550, 564, 565, 592, 594, 600, 615, 618, 654, 655, 666, 668, 676, 678, 679, 687, 723, 724, 736
    ]
    
    # Calculate the errors based on the first resolution
    errors_dict = {}
    base_results = results_dict[RESOLUTIONS[0]["resolution"]]
    for resolution in RESOLUTIONS[1:]:

        # Initialise error information
        res_kw = resolution["resolution"]
        comp_results = results_dict[res_kw]
        errors_dict[res_kw] = {"stress": [], "orientation": []}
        
        # Iterate through parameters
        for param_kw in PARAM_KW_LIST:
            
            # Get results with different resolutions but same parameters
            base_result = base_results[param_kw]
            comp_result = comp_results[param_kw]
            
            # Calculate stress error
            stress_error = get_stress(
                stress_list_1 = base_result[STRESS_FIELD],
                stress_list_2 = comp_result[STRESS_FIELD],
                strain_list_1 = base_result[STRAIN_FIELD],
                strain_list_2 = comp_result[STRAIN_FIELD],
                eval_strains  = EVAL_STRAINS,
            )

            # Get common grains
            # base_grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in base_result.keys() if "_phi_1" in key]
            # comp_grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in comp_result.keys() if "_phi_1" in key]
            # common_grain_ids = list(filter(lambda g_id: g_id in base_grain_ids, comp_grain_ids))

            # Calculate orientation error
            geodesic_grid = get_geodesics(
                grain_ids     = common_grain_ids,
                data_dict_1   = base_result,
                data_dict_2   = comp_result,
                strain_list_1 = base_result[STRAIN_FIELD],
                strain_list_2 = comp_result[STRAIN_FIELD],
                eval_strains  = EVAL_STRAINS
            )
            geodesic_error = np.average([np.sqrt(np.average([g**2 for g in gg])) for gg in geodesic_grid])

            # Compile errors
            errors_dict[res_kw]["stress"].append(stress_error)
            errors_dict[res_kw]["orientation"].append(geodesic_error)

    # Plot the evaluation times first
    resolution_list = [res["resolution"] for res in RESOLUTIONS]
    plt.figure(figsize=(5,5), dpi=200)
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    plt.xlabel("Resolution (µm)", fontsize=LABEL_SIZE)
    plt.ylabel("Average evaluation time (h)", fontsize=LABEL_SIZE)
    plt.bar(resolution_list, [r["time"] for r in RESOLUTIONS], color=COLOUR, zorder=3, width=WIDTH, edgecolor="black")
    padding = INCREMENT-WIDTH/2
    plt.xlim(min(resolution_list)-padding, max(resolution_list)+padding)
    plt.yscale("log")
    plt.ylim(1e-1,1e2)
    plt.xticks(resolution_list, fontsize=TICK_SIZE)
    plt.yticks(fontsize=TICK_SIZE)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    save_plot("results/mt_times.png")
            
    # Prepare error plotting
    stress_error_grid      = [errors_dict[res_kw]["stress"]*10 for res_kw in errors_dict.keys()]
    orientation_error_grid = [errors_dict[res_kw]["orientation"]*10 for res_kw in errors_dict.keys()]
    resolution_list        = [res["resolution"] for res in RESOLUTIONS[1:]]

    # Plot stress errors
    stress_error_grid = [[se*100 for se in stress_error_list] for stress_error_list in stress_error_grid]
    plot_boxplots(resolution_list, stress_error_grid, COLOUR)
    plt.xlabel("Resolution (µm)", fontsize=LABEL_SIZE)
    plt.ylabel(r"$E_{\sigma}$"+" (%)", fontsize=LABEL_SIZE)
    plt.xlim(min(resolution_list)-2.5, max(resolution_list)+2.5)
    plt.ylim(0, 1.6)
    # plt.gca().ticklabel_format(axis="y", style="sci", scilimits=(-2,-2))
    plt.gca().yaxis.major.formatter._useMathText = True
    plt.gca().yaxis.get_offset_text().set_fontsize(TICK_SIZE)
    save_plot("results/mt_se.png")

    # Plot geodesic errors
    orientation_error_grid = [[oe*180/np.pi for oe in orientation_error_list] for orientation_error_list in orientation_error_grid]
    plot_boxplots(resolution_list, orientation_error_grid, COLOUR)
    plt.xlabel("Resolution (µm)", fontsize=LABEL_SIZE)
    plt.ylabel(r"$E_{\Phi}$"+" (°)", fontsize=LABEL_SIZE)
    plt.xlim(min(resolution_list)-2.5, max(resolution_list)+2.5)
    plt.ylim(0, 0.5)
    # plt.gca().ticklabel_format(axis="y", style="sci", scilimits=(-2,-2))
    plt.gca().yaxis.major.formatter._useMathText = True
    plt.gca().yaxis.get_offset_text().set_fontsize(TICK_SIZE)
    save_plot("results/mt_ge.png")

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
                           vert=True, widths=WIDTH, whiskerprops=dict(linewidth=1), capprops=dict(linewidth=1))
    
    # Apply additional formatting to the boxplots
    for i in range(len(y_list_list)):
        patch = boxplots["boxes"][i]
        patch.set_facecolor(colour)
        patch.set_edgecolor("black")
        patch.set_linewidth(1)
        median = boxplots["medians"][i]
        median.set(color="black", linewidth=1)

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
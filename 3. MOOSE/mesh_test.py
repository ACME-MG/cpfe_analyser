"""
 Title:         Convergence
 Description:   Analyses the convergence of different mesh resolutions
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import os, math, numpy as np
import matplotlib.pyplot as plt
from __common__.plotter import save_plot
from __common__.analyse import get_geodesics, get_stress
from scipy.interpolate import splev, splrep, splder
from scipy.spatial.transform import Rotation

# Constants
RESOLUTIONS = [
    {"resolution": 5,  "ebsd_id": "ebsd_1", "colour": "silver"},
    {"resolution": 10, "ebsd_id": "ebsd_2", "colour": "orange"},
    {"resolution": 15, "ebsd_id": "ebsd_3", "colour": "gold"},
    {"resolution": 20, "ebsd_id": "ebsd_4", "colour": "brown"},
    {"resolution": 25, "ebsd_id": "ebsd_2", "colour": "red"},
    {"resolution": 30, "ebsd_id": "ebsd_3", "colour": "magenta"},
    {"resolution": 35, "ebsd_id": "ebsd_2", "colour": "purple"},
    {"resolution": 40, "ebsd_id": "ebsd_4", "colour": "blue"},
    {"resolution": 45, "ebsd_id": "ebsd_3", "colour": "cyan"},
    # {"resolution": 50, "ebsd_id": "ebsd_4", "colour": "green"},
]
PARAM_KW_LIST  = ["p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7"]
GRAIN_MAP      = "data/res_grain_map.csv"
SIM_PATH       = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim/2024-09-26 (617_s3_converge_5pct_8p)"
STRAIN_FIELD   = "average_strain"
STRESS_FIELD   = "average_stress"
EVAL_STRAINS   = np.linspace(0, 0.05, 50)

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
    common_grain_ids = [29, 72, 77, 81, 87, 97, 101, 114, 132, 154, 167, 189, 203, 237, 264, 265, 279, 284, 288, 289, 302, 314, 317, 326, 328, 352, 365, 376, 380, 381, 392, 422, 427, 432, 438, 447, 453, 455, 460, 486, 490, 493, 509, 522, 525, 530, 535, 546, 550, 564, 565, 592, 594, 600, 615, 618, 654, 655, 666, 668, 676, 678, 679, 687, 723, 724, 736]
    
    # Calculate the errors based on the first resolution
    errors_dict = {}
    base_results = results_dict[RESOLUTIONS[0]["resolution"]]
    for resolution in RESOLUTIONS[1:]:

        # Initialise error information
        res_kw = resolution["resolution"]
        comp_results = results_dict[res_kw]
        errors_dict[res_kw] = {"stress": [], "orientation": [], "reduced": []}
        
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
            errors_dict[res_kw]["reduced"].append(geodesic_error/np.pi + stress_error)
            
    # Prepare error plotting
    stress_error_grid = [errors_dict[res_kw]["stress"] for res_kw in errors_dict.keys()]
    orientation_error_grid = [errors_dict[res_kw]["orientation"] for res_kw in errors_dict.keys()]
    reduced_error_grid = [errors_dict[res_kw]["reduced"] for res_kw in errors_dict.keys()]
    resolution_list = [res["resolution"] for res in RESOLUTIONS[1:]]
    
    # Plot stress errors
    plot_boxplots(resolution_list, stress_error_grid, (0.6, 0.8, 1.0))
    plt.xlabel("Resolution (µm)", fontsize=24, labelpad=16)
    plt.ylabel(r"$E_{\sigma}$", fontsize=24, labelpad=16)
    plt.xlim(max(resolution_list)+2.5, min(resolution_list)-2.5)
    plt.ylim(0, 0.016)
    plt.gca().ticklabel_format(axis="y", style="sci", scilimits=(-3,-3))
    plt.gca().yaxis.major.formatter._useMathText = True
    plt.gca().yaxis.get_offset_text().set_fontsize(18)
    save_plot("results/plot_mt_se.png")

    # Plot geodesic errors
    plot_boxplots(resolution_list, orientation_error_grid, (1.0, 0.6, 0.0))
    plt.xlabel("Resolution (µm)", fontsize=24, labelpad=16)
    plt.ylabel(r"Average $E_{\Phi}$", fontsize=24, labelpad=16)
    plt.xlim(max(resolution_list)+2.5, min(resolution_list)-2.5)
    plt.ylim(0, 0.008)
    plt.gca().ticklabel_format(axis="y", style="sci", scilimits=(-3,-3))
    plt.gca().yaxis.major.formatter._useMathText = True
    plt.gca().yaxis.get_offset_text().set_fontsize(18)
    save_plot("results/plot_mt_ge.png")

    # Plot reduced errors
    plot_boxplots(resolution_list, reduced_error_grid, (0.6, 0.8, 1.0))
    plt.xlabel("Resolution (µm)", fontsize=24, labelpad=16)
    plt.ylabel(r"$E_{\Sigma}$", fontsize=24, labelpad=16)
    plt.xlim(max(resolution_list)+2.5, min(resolution_list)-2.5)
    plt.ylim(0, 0.018)
    plt.gca().ticklabel_format(axis="y", style="sci", scilimits=(-3,-3))
    plt.gca().yaxis.major.formatter._useMathText = True
    plt.gca().yaxis.get_offset_text().set_fontsize(18)
    save_plot("results/plot_mt_re.png")

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
                           vert=True, widths=4, whiskerprops=dict(linewidth=2), capprops=dict(linewidth=2))
    
    # Apply additional formatting to the boxplots
    for i in range(len(y_list_list)):
        patch = boxplots["boxes"][i]
        patch.set_facecolor(colour)
        patch.set_edgecolor("black")
        patch.set_linewidth(2)
        median = boxplots["medians"][i]
        median.set(color="black", linewidth=2)

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

def initialise_plot(x_label:str="", y_label:str="", x_max:float=None, y_max:float=None) -> None:
    """
    Initialises a plot
    
    Parameters:
    * `x_max`:   Upper limit for the x-axis
    * `y_max`:   Upper limit for the y-axis
    * `x_label`: Label for the x-axis
    * `y_label`: Label for the y-axis
    """
    plt.figure(figsize=(5,5))
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":")
    plt.xlabel(x_label, fontsize=11) if x_label != "" else None
    plt.ylabel(y_label, fontsize=11) if y_label != "" else None
    plt.xlim(0,x_max) if x_max != None else None
    plt.ylim(0,y_max) if y_max != None else None
    
def format_and_save_plot(plot_path:str, add_legend:bool=True, settings:dict={}) -> None:
    """
    Formats and saves a plot
    
    Parameters:
    * `plot_path`:  Path to save the plot
    * `add_legend`: Whether to add a legend
    * `settings`:   Settings for the legend
    """
    if add_legend:
        plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white", **settings)
    plt.savefig(plot_path)
    plt.cla()
    plt.clf()
    plt.close()

def euler_to_quat(euler:list) -> list:
    """
    Converts a set of euler-bunge angles (rads) into a quaternion

    Parameters:
    `euler`: The euler angle (rads)

    Returns the quaternion as a list (x,y,z,w)
    """
    euler_array = np.array(euler)
    rotation = Rotation.from_euler("zxz", euler_array, degrees=False)
    quat = rotation.as_quat()
    return list(quat)

def intervaluate(data_dict:dict, key:str) -> list:
    """
    Interpolates a list of values based on the strain and evaluates
    the interpolator using the defined strain list
    
    Parameters:
    * `data_dict`: Dictionary of values
    * `key`:       Key in dictionary to intervaluate
    
    Returns the evaluation from the interpolator as a list
    """
    interpolator = Interpolator(data_dict[STRAIN_FIELD], data_dict[key])
    evaluation = interpolator.evaluate(STRAIN_LIST)
    return evaluation

# The Interpolator Class
class Interpolator:

    def __init__(self, x_list:list, y_list:list, resolution:int=50, smooth:bool=False):
        """
        Class for interpolating two lists of values

        Parameters:
        * `x_list`:     List of x values
        * `y_list`:     List of y values
        * `resolution`: The resolution used for the interpolation
        * `smooth`:     Whether to smooth the interpolation
        """
        x_list, indices = np.unique(np.array(x_list), return_index=True)
        y_list = np.array(y_list)[indices]
        if len(x_list) > resolution:
            x_list = get_thinned_list(list(x_list), resolution)
            y_list = get_thinned_list(list(y_list), resolution)
        smooth_amount = resolution if smooth else 0
        self.spl = splrep(x_list, y_list, s=smooth_amount)
    
    def differentiate(self) -> None:
        """
        Differentiate the interpolator
        """
        self.spl = splder(self.spl)

    def evaluate(self, x_list:list) -> list:
        """
        Run the interpolator for specific values

        Parameters
        * `x_list`: The list of x values

        Returns the evaluated values
        """
        return list(splev(x_list, self.spl))

def get_thinned_list(unthinned_list:list, density:int) -> list:
    """
    Gets a thinned list

    Parameters:
    * `unthinned_list`: The list before thinning
    * `density`:        The goal density of the thinned list

    Returns the thinned list
    """
    src_data_size = len(unthinned_list)
    step_size = src_data_size / density
    thin_indexes = [math.floor(step_size*i) for i in range(1, density - 1)]
    thin_indexes = [0] + thin_indexes + [src_data_size - 1]
    thinned_list = [unthinned_list[i] for i in thin_indexes]
    return thinned_list

def csv_to_dict(csv_path:str, delimeter:str=",") -> dict:
    """
    Converts a CSV file into a dictionary
    
    Parameters:
    * `csv_path`:  The path to the CSV file
    * `delimeter`: The separating character
    
    Returns the dictionary
    """

    # Read all data from CSV (assume that file is not too big)
    csv_fh = open(csv_path, "r", encoding="utf-8-sig")
    csv_lines = csv_fh.readlines()
    csv_fh.close()

    # Initialisation for conversion
    csv_dict = {}
    headers = csv_lines[0].replace("\n", "").split(delimeter)
    csv_lines = csv_lines[1:]
    for header in headers:
        csv_dict[header] = []

    # Start conversion to dict
    for csv_line in csv_lines:
        csv_line_list = csv_line.replace("\n", "").split(delimeter)
        for i in range(len(headers)):
            value = csv_line_list[i]
            if value == "":
                continue
            try:
                value = float(value)
            except:
                pass
            csv_dict[headers[i]].append(value)
    
    # Convert single item lists to items and things multi-item lists
    for header in headers:
        if len(csv_dict[header]) == 1:
            csv_dict[header] = csv_dict[header][0]
    
    # Return
    return csv_dict
    
# Main function
if __name__ == "__main__":
    main()
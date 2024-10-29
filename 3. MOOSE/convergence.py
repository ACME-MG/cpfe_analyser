"""
 Title:         Convergence
 Description:   Analyses the convergence of different mesh resolutions
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += ["/home/janzen/code/moose_sim"]
import os, math, numpy as np
import matplotlib.pyplot as plt
from neml.math import rotations
from neml.cp import crystallography
from scipy.interpolate import splev, splrep, splder
from scipy.spatial.transform import Rotation
from moose_sim.helper.general import transpose
from moose_sim.analyse.pole_figure import IPF, get_lattice

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
    {"resolution": 50, "ebsd_id": "ebsd_4", "colour": "green"},
]
PARAM_KW_LIST = ["p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7"]
GRAIN_MAP     = "data/res_grain_map.csv"
SIM_PATH      = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim/2024-09-26 (617_s3_converge_5pct_8p)"
STRAIN_KEY    = "average_strain"
STRAIN_LIST   = [0.01, 0.02, 0.03, 0.04, 0.05]

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
    
    # Plot raw stress-strain data
    initialise_plot("Strain (mm/mm)", "Stress (MPa)", 0.05, 600)
    for i, resolution in enumerate(RESOLUTIONS):
        colour = resolution["colour"]
        res_kw = resolution["resolution"]
        settings = {"linewidth": 3} if i == 0 else {"linestyle": "dashed"}
        for param_kw in PARAM_KW_LIST:
            sum_dict = results_dict[res_kw][param_kw]
            plt.plot(sum_dict["average_strain"], sum_dict["average_stress"], colour, **settings)
        plt.plot([], [], label=f"{res_kw}µm", color=colour, **settings)
    format_and_save_plot("results/plot_ss.png")
    
    # Identify common grains across all meshes
    grain_ids_list = []
    for resolution in RESOLUTIONS:
        res_kw = resolution["resolution"]
        for param_kw in PARAM_KW_LIST:
            sum_dict = results_dict[res_kw][param_kw]
            grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in sum_dict.keys() if "_phi_1" in key]
            grain_ids_list.append(grain_ids)
    common_grain_ids = grain_ids_list[0]
    for grain_ids in grain_ids_list[1:]:
        common_grain_ids = list(filter(lambda g_id: g_id in grain_ids, common_grain_ids))
    
    # Plot certain common trajectories
    grain_ids = common_grain_ids[:5]
    get_trajectories = lambda dict : [transpose([dict[f"g{grain_id}_{phi}"][1:] for phi in ["phi_1", "Phi", "phi_2"]]) for grain_id in grain_ids]
    ipf = IPF(get_lattice("fcc"))
    direction = [1,0,0]
    for i, resolution in enumerate(RESOLUTIONS):
        res_kw = resolution["resolution"]
        colour = resolution["colour"]
        for param_kw in PARAM_KW_LIST:
            sum_dict = results_dict[res_kw][param_kw]
            trajectories = get_trajectories(sum_dict)
            if i == 0:
                ipf.plot_ipf_trajectory(trajectories, direction, "plot", {"color": colour, "linewidth": 2})
                ipf.plot_ipf_trajectory([[t[0]] for t in trajectories], direction, "scatter", {"color": colour, "s": 4**2})
                for exp_trajectory, grain_id in zip(trajectories, common_grain_ids):
                    ipf.plot_ipf_trajectory([[exp_trajectory[0]]], direction, "text", {"color": "black", "fontsize": 8, "s": grain_id})
            else:
                ipf.plot_ipf_trajectory(trajectories, direction, "plot", {"color": colour, "linewidth": 1, "zorder": 3})
                ipf.plot_ipf_trajectory([[t[0]] for t in trajectories], direction, "scatter", {"color": colour, "s": 3**2, "zorder": 3})
        settings = {"linewidth": 3} if i == 0 else {}
        plt.plot([], [], label=f"{res_kw}µm", color=colour, **settings)
    format_and_save_plot("results/plot_rt.png")
    
    # Calculate the errors based on the first resolution
    errors_dict = {}
    base_results = results_dict[RESOLUTIONS[0]["resolution"]]
    for resolution in RESOLUTIONS[1:]:
        res_kw = resolution["resolution"]
        comp_results = results_dict[res_kw]
        errors_dict[res_kw] = {"stress": [], "orientation": []}
        for param_kw in PARAM_KW_LIST:
            
            # Get results with different resolutions but same parameters
            base_result = base_results[param_kw]
            comp_result = comp_results[param_kw]
            
            # Calculate stress error
            base_stress_list = intervaluate(base_result, "average_grain_stress")
            comp_stress_list = intervaluate(comp_result, "average_grain_stress")
            stress_error = [abs((base_stress-comp_stress)) 
                            for base_stress, comp_stress in zip(base_stress_list, comp_stress_list)]
            stress_error = np.average(stress_error)*100/np.average(base_stress_list)

            # # Get common grains
            # base_grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in base_result.keys() if "_phi_1" in key]
            # comp_grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in comp_result.keys() if "_phi_1" in key]
            # common_grain_ids = list(filter(lambda g_id: g_id in base_grain_ids, comp_grain_ids))

            # Calculate orientation error for each common grain
            orientation_error = []
            for grain_id in common_grain_ids:
                base_phi_1_list = intervaluate(base_result, f"g{grain_id}_phi_1")
                base_Phi_list   = intervaluate(base_result, f"g{grain_id}_Phi")
                base_phi_2_list = intervaluate(base_result, f"g{grain_id}_phi_2")
                comp_phi_1_list = intervaluate(comp_result, f"g{grain_id}_phi_1")
                comp_Phi_list   = intervaluate(comp_result, f"g{grain_id}_Phi")
                comp_phi_2_list = intervaluate(comp_result, f"g{grain_id}_phi_2")
                for i in range(len(STRAIN_LIST)):
                    base_euler = [base_phi_1_list[i], base_Phi_list[i], base_phi_2_list[i]]
                    comp_euler = [comp_phi_1_list[i], comp_Phi_list[i], comp_phi_2_list[i]]
                    # orientation_error.append(get_cubic_misorientation(base_euler, comp_euler))
                    orientation_error.append(get_geodesic(euler_to_quat(base_euler), euler_to_quat(comp_euler)))
            orientation_error = np.average(orientation_error)
        
            # Compile errors
            errors_dict[res_kw]["stress"].append(stress_error)
            errors_dict[res_kw]["orientation"].append(orientation_error)
    
    # Prepare error plotting
    resolution_list = [res["resolution"] for res in RESOLUTIONS[1:]]
    
    # Plot stress-strain errors
    initialise_plot("Resolution (µm)", "Relative Error (%)")
    for res in resolution_list:
        plt.scatter([res]*len(errors_dict[res]["stress"]), errors_dict[res]["stress"], marker="o", s=6**2, alpha=0.50)
    average_errors = [np.average(errors_dict[res]["stress"]) for res in resolution_list]
    plt.plot(resolution_list, average_errors, color="black", linestyle="dashed", label="Average")
    plt.xlim(min(resolution_list)-2.5, max(resolution_list)+2.5)
    plt.ylim(0, 1.6)
    plt.gca().set_xticks(resolution_list)
    plt.gca().set_xticklabels(resolution_list)
    format_and_save_plot("results/err_ss.png", settings={"loc": "upper left"})

    # Plot orientation errors
    initialise_plot("Resolution (µm)", "Geodesic Distance (rads)")
    # initialise_plot("Resolution (µm)", "Misorientation (rads)")
    for res in resolution_list:
        plt.scatter([res]*len(errors_dict[res]["orientation"]), errors_dict[res]["orientation"], marker="o", s=6**2, alpha=0.50)
    average_errors = [np.average(errors_dict[res]["orientation"]) for res in resolution_list]
    plt.plot(resolution_list, average_errors, color="black", linestyle="dashed", label="Average")
    plt.xlim(min(resolution_list)-2.5, max(resolution_list)+2.5)
    # plt.ylim(0.002, 0.008)
    plt.ylim(0.0, 0.005)
    plt.gca().set_xticks(resolution_list)
    plt.gca().set_xticklabels(resolution_list)
    format_and_save_plot("results/err_rt.png", settings={"loc": "upper left"})

def get_cubic_misorientation(euler_1:list, euler_2:list):
    """
    Determines the misorientation of two sets of euler angles (rads);
    assumes cubic structure

    Parameters:
    * `euler_1`: The first euler angle
    * `euler_2`: The second euler angle
    
    Returns the misorientation angle
    """
    euler_1 = rotations.CrystalOrientation(*euler_1, angle_type="radians", convention="bunge")
    euler_2 = rotations.CrystalOrientation(*euler_2, angle_type="radians", convention="bunge")
    sym_group = crystallography.SymmetryGroup("432")
    mori = sym_group.misorientation(euler_1, euler_2)
    _, mori_angle = mori.to_axis_angle()
    return mori_angle

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

def get_geodesic(quat_1:list, quat_2:list) -> float:
    """
    Gets the geodesic distance between two quaternions angles
    
    Parameters:
    * `quat_1`: The first quaternion
    * `quat_2`: The second quaternion
    
    Returns the geodesic distance
    """
    quat_1 = np.array(quat_1)
    quat_2 = np.array(quat_2)
    quat_1 = quat_1 / np.linalg.norm(quat_1)
    quat_2 = quat_2 / np.linalg.norm(quat_2)
    dot_product = np.dot(quat_1, quat_2)
    dot_product = np.clip(dot_product, -1.0, 1.0)
    distance = np.arccos(np.abs(dot_product))
    return distance

def intervaluate(data_dict:dict, key:str) -> list:
    """
    Interpolates a list of values based on the strain and evaluates
    the interpolator using the defined strain list
    
    Parameters:
    * `data_dict`: Dictionary of values
    * `key`:       Key in dictionary to intervaluate
    
    Returns the evaluation from the interpolator as a list
    """
    interpolator = Interpolator(data_dict[STRAIN_KEY], data_dict[key])
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
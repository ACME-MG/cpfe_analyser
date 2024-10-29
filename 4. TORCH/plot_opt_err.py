"""
 Title:         Optimisation Test
 Description:   Tests the surrogate
 Author:        Janzen Choi

"""

# Libraries
import torch, math, numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import sys; sys.path += ["..", "/home/janzen/code/mms"]
from __common__.io import csv_to_dict
from __common__.general import transpose
from __common__.plotter import save_plot
from __common__.interpolator import intervaluate
from __common__.orientation import euler_to_quat, get_geodesic

# Paths
DIRECTORY = "617_s3_40um_lh6"
SUR_PATH = f"data/{DIRECTORY}/sm.pt"
MAP_PATH = f"data/{DIRECTORY}/map.csv"
EXP_PATH = "data/617_s3_exp.csv"
SIM_PATH = "data/summary.csv"

# Constants
EVAL_STRAINS = [0.05, 0.10, 0.15, 0.20, 0.25]
MAX_STRAIN = 0.1
PLOT_SIM = True

# Define model parameters
PARAM_LIST = [267.79, 122.95, 181.97, 328.73, 295.15, 238.77, 88.089, 7.65]
CAL_GRAIN_IDS = [243, 53, 240, 279, 256]
VAL_GRAIN_IDS = [54, 137, 271, 239, 207]
# VAL_GRAIN_IDS = [54, 137, 271, 239, 207, 68, 224, 265, 32, 56, 21, 115, 306, 101, 312]

# Main function
def main():

    # Get all results
    model = Model(SUR_PATH, MAP_PATH, EXP_PATH)
    res_dict = csv_to_dict(SIM_PATH) if PLOT_SIM else model.get_response(PARAM_LIST)
    exp_dict = csv_to_dict(EXP_PATH)

    # Get strain data
    strain_field = "average_strain" if PLOT_SIM else "strain_intervals"
    res_strain_list = res_dict[strain_field]
    exp_strain_list = exp_dict["strain_intervals"]

    # Get geodesic distances
    all_grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in res_dict.keys() if "_phi_1" in key]
    geodesic_grid = get_geodesics(CAL_GRAIN_IDS, res_dict, exp_dict, res_strain_list, exp_strain_list)
    remove_outlier = lambda data : [x for x in data if (q1 := np.percentile(data, 25)) - 1.5 *
                                    (iqr := np.percentile(data, 75) - q1) <= x <= q1 + 1.5 * iqr]
    print(geodesic_grid)
    geodesic_grid = [remove_outlier(geodesic_list) for geodesic_list in geodesic_grid]
    print(geodesic_grid)

    # Plot geodesic errors
    plt.figure(figsize=(6, 6))
    plot_boxplots(transpose(geodesic_grid), ["green"]*5)
    save_plot("results/boxplot.png")

def get_geodesics(grain_ids:list, data_dict_1:dict, data_dict_2:dict,
                  strain_list_1:list, strain_list_2:list) -> list:
    """
    Gets the geodesic distances from two sets of data
    
    Parameters:
    * `grain_ids`:    List of grain IDs to conduct evaluation
    * `data_dict_1`:  First set of data
    * `data_dict_2`:  Second set of data
    * `strain_list_1`: First list of strain values
    * `strain_list_2`: Second list of strain values
    
    Returns list of lists of geodesic distances
    """

    # Initialise
    geodesic_grid = []

    # Iterate through grain IDs
    for grain_id in grain_ids:

        # Get list of euler angles at specific strains
        quick_ie = lambda data_dict, strain_list : intervaluate_eulers(*[data_dict[f"g{grain_id}_{phi}"]
                   for phi in ["phi_1", "Phi", "phi_2"]], strain_list, EVAL_STRAINS)
        euler_list_1 = quick_ie(data_dict_1, strain_list_1)
        euler_list_2 = quick_ie(data_dict_2, strain_list_2)
        
        # Calculate geodesic distances of orientations at the same strains
        geodesic_list = []
        for euler_1, euler_2 in zip(euler_list_1, euler_list_2):
            quat_1 = euler_to_quat(euler_1)
            quat_2 = euler_to_quat(euler_2)
            geodesic = get_geodesic(quat_1, quat_2)
            geodesic_list.append(geodesic)
        geodesic_grid.append(geodesic_list)

    # Return list of lists of geodesic distances
    return geodesic_grid

def intervaluate_eulers(phi_1_list:list, Phi_list:list, phi_2_list:list,
                       strain_list:list, eval_strains:list) -> list:
    """
    Interpolates the euler-bunge (rads) components and evaluates
    the components and certain strain values

    Parameters:
    * `phi_1_list`:   List of phi_1 components
    * `Phi_list`:     List of Phi components
    * `phi_2_list`:   List of phi_2 components
    * `strain_list`:  List of strain values to interpolate
    * `eval_strains`: List of strain values to evaluate
    
    Returns the evaluated orientations as a list of euler-bunge values (in lists)
    """
    phi_1_list = intervaluate(strain_list, phi_1_list, eval_strains)
    Phi_list = intervaluate(strain_list, Phi_list, eval_strains)
    phi_2_list = intervaluate(strain_list, phi_2_list, eval_strains)
    euler_list = transpose([phi_1_list, Phi_list, phi_2_list])
    return euler_list

def linear(value:float, map:dict, mapper, index:int) -> float:
    """
    Linearly maps or unmaps a value

    Parameters:
    * `value`:  The value to be mapped / unmapped
    * `map`:    The mapping information
    * `mapper`: The mapping function handler
    * `index`:  The index of the map
    """
    return mapper(
        value = value,
        in_l  = map["in_l_bound"][index],
        in_u  = map["in_u_bound"][index],
        out_l = map["out_l_bound"][index],
        out_u = map["out_u_bound"][index],
    )

def linear_map(value:float, in_l:float, in_u:float, out_l:float, out_u:float) -> float:
    """
    Linearly maps a value

    Parameters:
    * `value`:  The value to be mapped
    * `in_l`:   The lower bound of the input
    * `in_u`:   The upper bound of the input
    * `out_l`:  The lower bound of the output
    * `out_u`:  The upper bound of the output

    Returns the mapped value
    """
    if in_l == in_u or out_l == out_u:
        return value
    factor = (out_u - out_l) / (in_u - in_l)
    return (value - in_l) * factor + out_l

def linear_unmap(value:float, in_l:float, in_u:float, out_l:float, out_u:float) -> float:
    """
    Linearly unmaps a value

    Parameters:
    * `value`:  The value to be unmapped
    * `in_l`:   The lower bound of the input
    * `in_u`:   The upper bound of the input
    * `out_l`:  The lower bound of the output
    * `out_u`:  The upper bound of the output

    Returns the unmapped value
    """
    if in_l == in_u or out_l == out_u:
        return value
    factor = (out_u - out_l) / (in_u - in_l)
    return (value - out_l) / factor + in_l

def get_sm_info(sm_path:str, map_path:str) -> tuple:
    """
    Loads the model and maps given two paths

    Parameters:
    * `sm_path`:    Path to the surrogate model
    * `map_path`:   Path to the map

    Returns the surrogate model, the input map, and the output map
    """

    # Load surrogate model
    model = torch.load(sm_path)
    model.eval()
    
    # Load maps
    input_map_dict, output_map_dict = {}, {}
    map_dict = csv_to_dict(map_path)
    num_inputs = map_dict["param_type"].count("input")
    for key in map_dict.keys():
        input_map_dict[key] = map_dict[key][:num_inputs]
        output_map_dict[key] = map_dict[key][num_inputs:]

    # Return everything
    return model, input_map_dict, output_map_dict

# Class for the model
class Model:
    
    def __init__(self, sm_path:str, map_path:str, exp_path:str) -> None:
        """
        Constructor for the surrogate model

        Parameters:
        * `sm_path`:  The path to the surrogate model
        * `map_path`: The path to the mapping methods
        * `exp_path`: The path to the experimental data
        """
        
        # Get model and maps
        self.model, self.input_map, self.output_map = get_sm_info(sm_path, map_path)
        
        # Extract experimental information
        exp_dict = csv_to_dict(exp_path)
        self.response_dict = {"strain": np.linspace(0, MAX_STRAIN, 101)}
        # self.response_dict = {"strain": exp_dict["strain_intervals"]}
        self.response_dict["strain_intervals"] = self.response_dict["strain"]
        for param_name in self.output_map["param_name"]:
            is_orientation = "phi_1" in param_name or "Phi" in param_name or "phi_2" in param_name
            initial_value = exp_dict[param_name][0] if is_orientation else 0.0
            self.response_dict[param_name] = [initial_value]

    def get_response(self, param_list:list) -> dict:
        """
        Gets the response of the model from the parameters

        Parameters:
        * `param_list`: List of parameters
        
        Returns the response as a dictionary
        """

        # Initialise
        response_dict = deepcopy(self.response_dict)

        # Get outputs and combine
        for strain in response_dict["strain"][1:]:
            output_dict = self.get_output(param_list + [strain])
            if output_dict == None:
                return
            for key in output_dict.keys():
                response_dict[key].append(output_dict[key])
        
        # Adjust and return
        if "average_stress" in response_dict.keys():
            response_dict["stress"] = response_dict.pop("average_stress")
        return response_dict

    def get_output(self, input_list:list) -> dict:
        """
        Gets the outputs of the surrogate model

        Parameters:
        * `input_list`: The list of raw input values

        Returns the outputs
        """
        
        # Process inputs
        processed_input_list = []
        for i in range(len(input_list)):
            try:
                input_value = math.log(input_list[i]) / math.log(self.input_map["base"][i])
                input_value = linear(input_value, self.input_map, linear_map, i)
            except ValueError:
                return None
            processed_input_list.append(input_value)
        
        # Get raw outputs and process
        input_tensor = torch.tensor(processed_input_list)
        with torch.no_grad():
            output_list = self.model(input_tensor).tolist()
        for i in range(len(output_list)):
            try:
                output_list[i] = linear(output_list[i], self.output_map, linear_unmap, i)
                output_list[i] = math.pow(self.output_map["base"][i], output_list[i])
            except ValueError:
                return None
        
        # Return the dictionary of outputs
        output_dict = dict(zip(self.output_map["param_name"], output_list))
        return output_dict

def plot_boxplots(y_list_list:list, colours:list) -> None:
    """
    Plots several boxplots together

    Parameters:
    * `y_list_list`: List of data lists
    * `colours`:     List of colours
    """

    # Plot boxplots
    boxplots = plt.boxplot(y_list_list, patch_artist=True, vert=True, widths=0.7)
    
    # Add scattered data and format each boxplot
    for i, (y_list, colour) in enumerate(zip(y_list_list, colours)):

        # Format boxplot face
        patch = boxplots["boxes"][i]
        patch.set_facecolor(colour)
        patch.set_alpha(0.4)
        patch.set_edgecolor("black")
        patch.set_linewidth(1.5)

        # Format median line
        median = boxplots["medians"][i]
        median.set(color="black", linewidth=1)

        # Add scattered data
        x_list = np.random.normal(i+1, 0.04, size=len(y_list))
        plt.scatter(x_list, y_list, s=4**2, color=colour, edgecolors="black")

# Calls the main function
if __name__ == "__main__":
    main()

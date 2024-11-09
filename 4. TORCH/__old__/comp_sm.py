"""
 Title:         Comparative Plotter for the surrogate model's errors
 Description:   Compares the errors of multiple surrogate models using a plot
 Author:        Janzen Choi

"""

# Libraries
import torch, math, numpy as np, os
import matplotlib.pyplot as plt
from copy import deepcopy
import sys; sys.path += ["..", "/home/janzen/code/mms"]
from __common__.io import csv_to_dict
from __common__.general import transpose, flatten
from __common__.plotter import save_plot
from __common__.interpolator import intervaluate
from __common__.orientation import euler_to_quat, get_geodesic

# Surrogate paths
SM_PATH  = "data/test_intervals"
EXP_PATH = "data/test_intervals/exp.csv"
SM_FILE  = "sm.pt"
MAP_FILE = "map.csv"

# Validation paths
VAL_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim/2024-11-04 (617_lh2_40um_sm)"
# VAL_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim/2024-10-08 (617_s3_lh_40um_sm)"
PRM_FILE = "params.txt"

# Other constants
PARAM_NAMES  = [f"cp_lh_{i}" for i in range(2)] + ["cp_tau_0", "cp_n"]
MAX_STRAIN   = 0.1
EVAL_STRAINS = [0.01*(i+1) for i in range(10)]
GRAIN_IDS    = [243, 53, 240, 279, 256]

def main():

    # Get surrogate paths
    interval_list = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48]
    sm_path_list = [f"{SM_PATH}/617_s3_lh2_{i}i" for i in interval_list]

    # Get validation paths
    val_path_list = [f"{VAL_PATH}/{path}" for path in os.listdir(VAL_PATH)]
    param_dict_list = [read_params(f"{path}/{PRM_FILE}") for path in val_path_list]

    # Initialise average error lists
    average_stress_error_list = []
    average_geodesic_error_list = []

    # Iterate through surrogate models
    for sm_path in sm_path_list:

        # Get model response
        model = Model(
            sm_path  = f"{sm_path}/{SM_FILE}",
            map_path = f"{sm_path}/{MAP_FILE}",
            exp_path = EXP_PATH
        )

        # Initialise error lists
        stress_error_list = []
        geodesic_error_list = []

        # Iterate through validation parameters
        for val_path, param_dict in zip(val_path_list, param_dict_list):

            # Get simulated response
            sim_dict = csv_to_dict(f"{val_path}/summary.csv")
            sim_strain_list = sim_dict["average_strain"]
            sim_stress_list = sim_dict["average_stress"]

            # Get surrogate model response
            param_values = [param_dict[param_name] for param_name in PARAM_NAMES]
            res_dict = model.get_response(param_values)
            res_strain_list = res_dict["strain"]
            res_stress_list = res_dict["stress"]

            # Calculate stress error
            sim_stress = intervaluate(sim_strain_list, sim_stress_list, EVAL_STRAINS)
            res_stress = intervaluate(res_strain_list, res_stress_list, EVAL_STRAINS)
            stress_error = np.average([abs((ss-rs)/ss) for ss, rs in zip(sim_stress, res_stress)])
            stress_error_list.append(stress_error)

            # Calculate geodesic error
            geodesic_grid = get_geodesics(GRAIN_IDS, res_dict, sim_dict, res_strain_list, sim_strain_list)
            geodesic_error = np.average(flatten(geodesic_grid))
            geodesic_error_list.append(geodesic_error)

        # Average the errors
        average_stress_error_list.append(np.average(stress_error_list)*100)
        average_geodesic_error_list.append(np.average(geodesic_error_list))
    
    # Plot stress errors
    plt.figure(figsize=(5, 5))
    plt.grid(True)
    plt.plot(interval_list, average_stress_error_list, marker="o")
    plt.xlabel("Number of Evaluations")
    plt.ylabel("Stress Relative Error (%)")
    plt.xlim(0, 50)
    plt.ylim(0, 6)
    save_plot("results/stress_error.png")

    # Plot geodesic errors
    plt.figure(figsize=(5, 5))
    plt.grid(True)
    plt.plot(interval_list, average_geodesic_error_list, marker="o")
    plt.xlabel("Number of Evaluations")
    plt.ylabel("Geodesic Error (rads)")
    plt.xlim(0, 50)
    plt.ylim(0, 0.0010)
    plt.gca().ticklabel_format(axis="y", style="sci", scilimits=(-4,-4))
    plt.gca().yaxis.major.formatter._useMathText = True
    save_plot("results/geodesic_error.png")

def get_geodesics(grain_ids:list, data_dict_1:dict, data_dict_2:dict,
                  strain_list_1:list, strain_list_2:list) -> list:
    """
    Gets the geodesic distances from two sets of data
    
    Parameters:
    * `grain_ids`:     List of grain IDs to conduct evaluation
    * `data_dict_1`:   First set of data
    * `data_dict_2`:   Second set of data
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

# Calls the main function
if __name__ == "__main__":
    main()


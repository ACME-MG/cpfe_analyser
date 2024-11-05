"""
 Title:         Optimisation Test
 Description:   Tests the surrogate
 Author:        Janzen Choi

"""

# Libraries
import torch, math, numpy as np
from copy import deepcopy
import sys; sys.path += ["..", "/home/janzen/code/mms"]
from __common__.io import csv_to_dict
from __common__.general import transpose
from __common__.pole_figure import get_lattice, IPF
from __common__.plotter import define_legend, save_plot, Plotter

# Paths
DIRECTORY = "617_s3_40um_lh2"
SUR_PATH = f"data/{DIRECTORY}/1_sm.pt"
MAP_PATH = f"data/{DIRECTORY}/1_map.csv"
EXP_PATH = f"data/{DIRECTORY}/exp.csv"
SIM_PATH = f"data/{DIRECTORY}/1_opt.csv"

# Constants
EVAL_STRAINS = [0.05, 0.10, 0.15, 0.20, 0.25]
MAX_STRAIN = 0.1
PLOT_SIM = True

# Define model parameters
PARAM_LIST = [396.17, 32.234, 79.898, 7.5364]
# CAL_GRAIN_IDS = [145, 243, 111, 49, 136]
CAL_GRAIN_IDS = [207, 79, 164, 167, 309]
VAL_GRAIN_IDS = []

# Main function
def main():

    # Get all results
    model = Model(SUR_PATH, MAP_PATH, EXP_PATH)
    res_dict = csv_to_dict(SIM_PATH) if PLOT_SIM else model.get_response(PARAM_LIST)
    exp_dict = csv_to_dict(EXP_PATH)
    get_trajectories = lambda dict, grain_ids : [transpose([dict[f"g{grain_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]) for grain_id in grain_ids]

    # Initialise IPF
    ipf = IPF(get_lattice("fcc"))
    direction = [1,0,0]

    # Plot experimental reorientation trajectories
    exp_trajectories = get_trajectories(exp_dict, CAL_GRAIN_IDS+VAL_GRAIN_IDS)
    ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": "silver", "linewidth": 2})
    ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": "silver", "head_width": 0.01, "head_length": 0.015})
    ipf.plot_ipf_trajectory([[et[0]] for et in exp_trajectories], direction, "scatter", {"color": "silver", "s": 8**2})
    for exp_trajectory, grain_id in zip(exp_trajectories, CAL_GRAIN_IDS+VAL_GRAIN_IDS):
        ipf.plot_ipf_trajectory([[exp_trajectory[0]]], direction, "text", {"color": "black", "fontsize": 8, "s": grain_id})

    # Plot calibration reorientation trajectories
    cal_trajectories = get_trajectories(res_dict, CAL_GRAIN_IDS)
    ipf.plot_ipf_trajectory(cal_trajectories, direction, "plot", {"color": "green", "linewidth": 1, "zorder": 3})
    ipf.plot_ipf_trajectory(cal_trajectories, direction, "arrow", {"color": "green", "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
    ipf.plot_ipf_trajectory([[ct[0]] for ct in cal_trajectories], direction, "scatter", {"color": "green", "s": 6**2, "zorder": 3})

    # Plot calibration reorientation trajectories
    val_trajectories = get_trajectories(res_dict, VAL_GRAIN_IDS)
    ipf.plot_ipf_trajectory(val_trajectories, direction, "plot", {"color": "red", "linewidth": 1, "zorder": 3})
    ipf.plot_ipf_trajectory(val_trajectories, direction, "arrow", {"color": "red", "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
    ipf.plot_ipf_trajectory([[vt[0]] for vt in val_trajectories], direction, "scatter", {"color": "red", "s": 6**2, "zorder": 3})

    # Save IPF
    define_legend(["silver", "green", "red"], ["Experimental", "Calibration", "Validation"], ["scatter", "line", "line"])
    save_plot("results/plot_rt.png")

    # Plot stress-strain curve
    if PLOT_SIM:
        res_dict["strain"] = res_dict["average_grain_strain"]
        res_dict["stress"] = res_dict["average_stress"]
    if "stress" in res_dict.keys():
        plotter = Plotter("strain", "stress", "mm/mm", "MPa")
        plotter.prep_plot()
        plotter.scat_plot(exp_dict, "silver", "Experimental")
        plotter.line_plot(res_dict, "green", "Calibration")
        plotter.set_legend()
        save_plot("results/plot_ss.png")

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

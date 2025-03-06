"""
 Title:         Plot Texture
 Description:   Plots the texture of the optimised simulations at certain strains
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import numpy as np
from __common__.general import transpose
from __common__.io import csv_to_dict
from __common__.plotter import save_plot
from __common__.pole_figure import get_lattice, PF
from __common__.interpolator import intervaluate

# Paths
EXP_PATH     = "data/617_s3_40um_exp.csv"
RESULTS_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/"
SIM_PATH     = f"{RESULTS_PATH}/asmbo/2025-02-28 (vh_pinned_sm8_i29)/250228104126_i28_simulate"

# Plotting parameters
EXP_COLOUR   = "silver"
SIM_COLOUR   = "tab:cyan" # VH
# SIM_COLOUR   = "tab:orange" # LH2
# SIM_COLOUR   = "tab:purple" # LH6
EVAL_STRAINS = [0.30]

def main():
    """
    Main function
    """

    # Initialise
    get_grain_ids    = lambda dict : [int(key.replace("g","").replace("_phi_1","")) for key in dict.keys() if "_phi_1" in key]
    get_trajectories = lambda dict, grain_ids : [transpose([dict[f"g{grain_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]) for grain_id in grain_ids]

    # Read experimental data
    exp_dict = csv_to_dict(EXP_PATH)
    exp_strain_list = exp_dict["strain_intervals"]
    exp_grain_ids = get_grain_ids(exp_dict)
    exp_trajectories = get_trajectories(exp_dict, exp_grain_ids)

    # Read simulation data
    sim_dict = csv_to_dict(f"{SIM_PATH}/summary.csv")
    sim_strain_list = sim_dict["average_strain"]
    sim_grain_ids = get_grain_ids(sim_dict)
    sim_trajectories = get_trajectories(sim_dict, sim_grain_ids)

    # Plot orientations for each evaluation strain
    for i, eval_strain in enumerate(EVAL_STRAINS):
        
        # Get orientations
        exp_orientations = [intervaluate_orientation(exp_strain_list, et, eval_strain) for et in exp_trajectories]
        sim_orientations = [intervaluate_orientation(sim_strain_list, st, eval_strain) for st in sim_trajectories]
        
        # Plot orientations for each direction
        # for direction in [[1,0,0], [1,1,0], [1,1,1]]:
        for direction in [[1,1,1]]:
            dir_str = "".join([str(d) for d in direction])
            
            # pf.plot_pf(exp_orientations, direction, colour=EXP_COLOUR)
            
            pf = PF(get_lattice("fcc"))
            pf.plot_pf(sim_orientations, direction, colour=SIM_COLOUR)
            save_plot(f"results/pf_{dir_str}_a{i+1}")
            
            pf = PF(get_lattice("fcc"))
            pf.plot_pf_density(sim_orientations, direction, colour=SIM_COLOUR)
            save_plot(f"results/pf_{dir_str}_b{i+1}")

def intervaluate_orientation(strain_list:list, orientation_list:list, strain:float) -> list:
    """
    Interpolates and evaluates from a list of orientations

    Parameters:
    * `strain_list`:      List of strain values
    * `orientation_list`: List of orientation components
    * `strain`:           Strain value to evaluate

    Returns the intervaluated orientation
    """
    components = transpose(orientation_list)
    new_components = [intervaluate(strain_list, component, strain) for component in components]
    return transpose(new_components)

# Main function
if __name__ == "__main__":
    main()

"""
 Title:         Plot Texture
 Description:   Plots the texture of the optimised simulations at certain strains
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import numpy as np
from __common__.general import transpose, round_sf
from __common__.io import csv_to_dict
from __common__.plotter import save_plot
from __common__.pole_figure import get_lattice, PFD, PF
from __common__.interpolator import intervaluate

# Paths
EXP_PATH     = "data/617_s3_40um_exp.csv"
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo"
MOOSE_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim"
RESULTS_PATH = "results/pf"

# Simulation information
SIM_INFO_LIST = [
    {"path": f"{MOOSE_PATH}/2025-03-15 (617_s3_vh_x_hr)", "colour": (1.0, 0.2, 0.2), "alpha": (0.2, 1.0)},
    # {"path": f"{ASMBO_PATH}/2025-03-03 (vh_pinned_sm8_i43)/250303022723_i26_simulate", "colour": (1.0, 0.2, 0.2), "alpha": (0.2, 1.0)},
    # {"path": f"{MOOSE_PATH}/2025-03-15 (617_s3_vh_c44d2)", "colour": (1.0, 0.2, 0.2), "alpha": (0.2, 1.0)},
    # {"path": f"{ASMBO_PATH}/2025-03-14 (vh_sbc_sm8_i37)/250314054502_i24_simulate", "colour": (1.0, 0.2, 0.2), "alpha": (0.2, 1.0)},
    # {"path": f"{ASMBO_PATH}/2025-03-03 (vh_pinned_sm8_i43)/250303022723_i26_simulate", "colour": (1.0, 0.2, 0.2), "alpha": (0.2, 1.0)},
    # {"path": f"{ASMBO_PATH}/2025-03-06 (vh_pin2_sm8_i40)/250306051551_i27_simulate", "colour": (1.0, 0.2, 0.2), "alpha": (0.6, 1.0)},
    # {"path": f"{ASMBO_PATH}/2025-03-09 (vh_pin2_sm8_i22)/250309175121_i10_simulate", "colour": (1.0, 0.2, 0.2), "alpha": (0.6, 1.0)}
]
EXP_COLOUR = {"colour": (0, 0, 0), "alpha": (0, 0.4)}

# Other parameters
EVAL_STRAINS = [0.30]
NUM_LEVELS   = 5

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
    for sim_info in SIM_INFO_LIST:
        sim_info["dict"] = csv_to_dict(f"{sim_info['path']}/summary.csv")
        sim_info["strain"] = sim_info["dict"]["average_strain"]
        sim_grain_ids = get_grain_ids(sim_info["dict"])
        sim_info["trajectories"] = get_trajectories(sim_info["dict"], sim_grain_ids)

    # Plot orientations for each evaluation strain
    for i, eval_strain in enumerate(EVAL_STRAINS):

        # Get experimental orientations
        exp_orientations = [intervaluate_orientation(exp_strain_list, et, eval_strain) for et in exp_trajectories]
        
        # Get simulated orientations
        sim_orientations_list = []
        for sim_info in SIM_INFO_LIST:
            sim_orientations = [intervaluate_orientation(sim_info["strain"], st, eval_strain) for st in sim_info["trajectories"]]
            sim_orientations_list.append(sim_orientations)
        
        # Plot orientations for each direction
        for direction in [[1,0,0], [1,1,0], [1,1,1]]:
            dir_str = "".join([str(d) for d in direction])
            
            # # Plot raw experimental texture
            # pf = PF(get_lattice("fcc"))
            # pf.plot_pf(exp_orientations, direction, colour=EXP_COLOUR["colour"])
            # save_plot(f"{RESULTS_PATH}/exp_{dir_str}_e{i+1}")

            # # Plot raw simulated texture
            # for j, sim_orientations in enumerate(sim_orientations_list):
            #     pf = PF(get_lattice("fcc"))
            #     pf.plot_pf(sim_orientations, direction, colour=sim_info["colour"])
            #     save_plot(f"{RESULTS_PATH}/sim{j+1}_{dir_str}_e{i+1}")

            # Plot contoured texture together
            pfd = PFD(get_lattice("fcc"))
            pfd.plot_pfd(exp_orientations, direction, levels=NUM_LEVELS, colour=EXP_COLOUR["colour"], alpha_limits=EXP_COLOUR["alpha"])
            for sim_orientations in sim_orientations_list:
                pfd.plot_pfd(sim_orientations, direction, levels=NUM_LEVELS, colour=sim_info["colour"], alpha_limits=sim_info["alpha"])
            save_plot(f"{RESULTS_PATH}/con_{dir_str}_e{i+1}")
        
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

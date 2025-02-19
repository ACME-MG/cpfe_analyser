"""
 Title:         Plot Optimised Simulations
 Description:   Plots the response of the optimised simulations
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import matplotlib.pyplot as plt
import numpy as np
from __common__.io import csv_to_dict
from __common__.general import transpose, round_sf
from __common__.pole_figure import get_lattice, IPF
from __common__.plotter import save_plot, Plotter
from __common__.analyse import get_geodesics, get_stress

# Experimental Information
EXP_PATH = "data/617_s3_40um_exp.csv"
EXP_COLOUR = "silver"
EXP_EBSD_ID = "ebsd_4"

# Simulation Information
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo"
MOOSE_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim"
SIM_INFO_LIST = [

    # Voce Hardening Model
    # {"label": "Low-fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-02-02 (vh_sm8_i72)/250202092030_i59_simulate"},
    # {"label": "High-fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/2025-02-04 (617_s3_10um_vh)"},
    
    # Two-coefficient latent hardening model
    # {"label": "Low-fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-02-06 (lh2_sm8_i15)/250206013451_i12_simulate"},
    # {"label": "High-fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/2025-02-09 (617_s3_10um_lh2)"},
    {"label": "Low-fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-02-19 (lh2_sm8_i17)/250218212104_i8_simulate"},
    # {"label": "High-fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/"},
    
    # Six-coefficient latent hardening model
    # {"label": "Low-fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-02-14 (lh6_sm16_i34)/250213134340_i26_simulate"},
    # {"label": "High-fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/2025-02-14 (617_s3_10um_lh6)"},
    
    # All low-fidelity
    # {"label": "VH",   "ebsd_id": "ebsd_4", "colour": "tab:cyan",   "path": f"{ASMBO_PATH}/2025-02-02 (vh_sm8_i72)/250202092030_i59_simulate"},
    # {"label": "LH2",  "ebsd_id": "ebsd_4", "colour": "tab:orange", "path": f"{ASMBO_PATH}/2025-02-06 (lh2_sm8_i15)/250206013451_i12_simulate"},
    # {"label": "LH6",  "ebsd_id": "ebsd_4", "colour": "tab:purple", "path": f"{ASMBO_PATH}/2025-02-14 (lh6_sm16_i34)/250213134340_i26_simulate"},
    
    # All high-fidelity
    {"label": "VH",   "ebsd_id": "ebsd_2", "colour": "tab:cyan",   "path": f"{MOOSE_PATH}/2025-02-04 (617_s3_10um_vh)"},
    {"label": "LH2",  "ebsd_id": "ebsd_2", "colour": "tab:orange", "path": f"{MOOSE_PATH}/2025-02-09 (617_s3_10um_lh2)"},
    {"label": "LH6",  "ebsd_id": "ebsd_2", "colour": "tab:purple", "path": f"{MOOSE_PATH}/2025-02-14 (617_s3_10um_lh6)"},
]
for si in SIM_INFO_LIST:
    si["data"] = csv_to_dict(f"{si['path']}/summary.csv")

# Grain IDs
GRAIN_IDS = [
    [51, 56, 72, 80, 126, 223, 237, 262], # Calibration
    [44, 50, 60, 178, 190, 207, 278, 299] # Validation
]

# Other Constants
STRAIN_FIELD = "average_strain"
STRESS_FIELD = "average_stress"
RES_DATA_MAP = "data/res_grain_map.csv"
# SPACING      = -2.25
SPACING      = -6.25

# Script parameters
SHOW_GRAIN_ID = False

# Main function
def main():

    # Prepare experimental data
    exp_dict = csv_to_dict(EXP_PATH)
    eval_strains = np.linspace(0, exp_dict["strain_intervals"][-1], 32)

    # Plot reorientation trajectories
    for i, grain_ids in enumerate(GRAIN_IDS):

        # Initialise IPF
        ipf = IPF(get_lattice("fcc"))
        direction = [1,0,0]
        get_trajectories = lambda dict, g_ids : [transpose([dict[f"g{g_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]) for g_id in g_ids]

        # Plot experimental reorientation trajectories
        exp_trajectories = get_trajectories(exp_dict, grain_ids)
        ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": EXP_COLOUR, "linewidth": 3})
        ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": EXP_COLOUR, "head_width": 0.01, "head_length": 0.015})
        ipf.plot_ipf_trajectory([[et[0]] for et in exp_trajectories], direction, "scatter", {"color": EXP_COLOUR, "s": 8**2})

        # Iterate through simulations
        for si in SIM_INFO_LIST:

            # Get simulated reorientation trajectories
            sim_grain_ids = [get_sim_grain_id(grain_id, si["ebsd_id"]) for grain_id in grain_ids]
            sim_trajectories = get_trajectories(si["data"], sim_grain_ids)
            if 44 in sim_grain_ids:
                sim_trajectories[sim_grain_ids.index(44)] = sim_trajectories[sim_grain_ids.index(44)][:-9] 
            
            # Plot simulated reorientation trajectories
            sim_colour = si["colour"]
            ipf.plot_ipf_trajectory(sim_trajectories, direction, "plot", {"color": sim_colour, "linewidth": 2, "zorder": 3})
            ipf.plot_ipf_trajectory(sim_trajectories, direction, "arrow", {"color": sim_colour, "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
            ipf.plot_ipf_trajectory([[st[0]] for st in sim_trajectories], direction, "scatter", {"color": sim_colour, "s": 6**2, "zorder": 3})

            # Plot grain IDs
            if SHOW_GRAIN_ID:
                for exp_trajectory, grain_id in zip(exp_trajectories, grain_ids):
                    ipf.plot_ipf_trajectory([[exp_trajectory[0]]], direction, "text", {"color": "blue", "fontsize": 8, "s": grain_id, "zorder": 3})

        # Add geodesic errors to legend
        ge_list = get_geodesic_errors(SIM_INFO_LIST, exp_dict, eval_strains, grain_ids)
        add_supp_legend(ge_list, SPACING)
        
        # Save IPF plot
        save_plot(f"results/plot_opts_rt_{i+1}.png")

    # Initialise stress-strain plot
    plotter = Plotter("strain", "stress", "mm/mm", "MPa")
    plotter.prep_plot(size=14)
    plotter.set_limits((0,0.5), (0,1600))

    # Plot stress-strain data
    plt.scatter(exp_dict["strain"], exp_dict["stress"], color=EXP_COLOUR, s=8**2)
    for si in SIM_INFO_LIST:
        plt.plot(si["data"][STRAIN_FIELD], si["data"][STRESS_FIELD], color=si["colour"], linewidth=3)

    # Add stress errors to legend
    sim_dict_list = [si["data"] for si in SIM_INFO_LIST]
    se_list = get_stress_errors(sim_dict_list, exp_dict, eval_strains)
    add_supp_legend(se_list, SPACING)

    # Format and save
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    save_plot("results/plot_opts_ss.png")

def add_supp_legend(error_list:list, spacing:float=-5.5) -> None:
    """
    Adds supplementary information to a legend that is aligned
    horizontally to the main keys; note that this immediately
    adds the legend to the current axis

    Parameters:
    * `error_list`: Error values to add as supplementary information
    * `spacing`:    Spacing between main keys and information
    """

    # Define main keys of the legend
    handles = [plt.scatter([], [], color=EXP_COLOUR, label="Experiment", s=8**2)]
    handles += [plt.plot([], [], color=si["colour"], label=si['label'], linewidth=3)[0] for si in SIM_INFO_LIST]

    # Define supplementary information
    se_label_list = [" "] + [f"({round_sf(error, 3)})" for error in error_list]
    handles += [plt.scatter([], [], color="white", label=se_label, marker="o", s=0) for se_label in se_label_list]

    # Add legend
    legend = plt.legend(handles=handles, ncol=2, columnspacing=spacing, framealpha=1, edgecolor="black",
                        fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    plt.gca().add_artist(legend)

def get_stress_errors(sim_dict_list:list, exp_dict:dict, eval_strains:list) -> list:
    """
    Calculates the stress errors of a list of simulations relative to experimental data

    Parameters:
    * `sim_dict_list`: The list of dictionaries of simulation results
    * `exp_dict`:      The dictionary of experimental data
    * `eval_strains`:  The strains to conduct the error evaluations
    
    Returns the stress errors as a list
    """
    stress_error_list = []
    for sim_dict in sim_dict_list:
        stress_error = get_stress(
            stress_list_1 = exp_dict["stress"],
            stress_list_2 = sim_dict["average_stress"],
            strain_list_1 = exp_dict["strain"],
            strain_list_2 = sim_dict["average_strain"],
            eval_strains  = eval_strains
        )
        stress_error_list.append(stress_error)
    return stress_error_list

def get_geodesic_errors(sim_info_list:list, exp_dict:dict, eval_strains:list, grain_ids:list) -> tuple:
    """
    Calculates the errors of a list of simulations relative to experimental data

    Parameters:
    * `sim_info_list`: The list of dictionaries of simulation results
    * `exp_dict`:      The dictionary of experimental data
    * `eval_strains`:  The strains to conduct the error evaluations
    * `grain_ids`:     The list of grain IDs
    
    Returns the geodesic errors as a list
    """

    # Iterate through simulations
    geodesic_error_list = []
    for si in sim_info_list:

        # Convert grain IDs
        eval_grain_ids = []
        sim_grain_ids = [get_sim_grain_id(grain_id, si["ebsd_id"]) for grain_id in grain_ids]
        sim_dict = {}
        for grain_id, sim_grain_id in zip(grain_ids, sim_grain_ids):
            if sim_grain_id == -1:
                continue
            for phi in ["phi_1", "Phi", "phi_2"]:
                sim_dict[f"g{grain_id}_{phi}"] = si["data"][f"g{sim_grain_id}_{phi}"]
            eval_grain_ids.append(grain_id)

        # Calculate geodesic errors
        geodesic_grid = get_geodesics(
            grain_ids     = eval_grain_ids,
            data_dict_1   = sim_dict,
            data_dict_2   = exp_dict,
            strain_list_1 = si["data"]["average_strain"],
            strain_list_2 = exp_dict["strain_intervals"],
            eval_strains  = eval_strains
        )

        # Compile geodesic errors
        geodesic_error = np.average([np.average(geodesic_list) for geodesic_list in geodesic_grid])
        geodesic_error_list.append(geodesic_error)
    
    # Return
    return geodesic_error_list

def get_sim_grain_id(exp_grain_id:int, ebsd_id:str) -> int:
    """
    Maps the experimental to simulated grain ID

    Parameters:
    * `exp_grain_id`: The grain ID from the experimental data
    * `ebsd_id`:      The origin of the EBSD map used to run the simulation

    Returns the corresponding grain ID in the simulated data
    """
    grain_map = csv_to_dict(RES_DATA_MAP)
    exp_ebsd_ids = grain_map[EXP_EBSD_ID]
    sim_ebsd_ids = grain_map[ebsd_id]
    return int(sim_ebsd_ids[exp_ebsd_ids.index(exp_grain_id)])

# Calls the main function
if __name__ == "__main__":
    main()

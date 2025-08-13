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
from __common__.general import transpose
from __common__.pole_figure import get_lattice, IPF
from __common__.plotter import save_plot, Plotter
from __common__.analyse import get_geodesics, get_stress
from sm_error import find_knee_point

# Experimental Information
EXP_PATH = "data/617_s3_40um_exp.csv"
EXP_COLOUR = "silver"
EXP_EBSD_ID = "ebsd_4"

# Simulation Information
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo"
MOOSE_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/moose_sim"
SIM_INFO_LIST = [
    
    # VH Model
    # {"label": "Run 1", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310145710_i22_simulate"},
    # # {"label": "Run 1", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-03-09 (vh_pin2_sm8_i25)/250308143546_i4_simulate"},
    # {"label": "Run 2", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310145710_i22_simulate"},
    # {"label": "Run 3", "ebsd_id": "ebsd_4", "colour": "tab:red",   "path": f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)/250318014435_i21_simulate"},
    # {"label": "Run 4", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-03-25 (vh_x_sm8_i31)/250325072901_i16_simulate"},
    # {"label": "Run 5", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310161708_i25_simulate"},
    {"label": "Surrogate", "ebsd_id": "ebsd_4", "colour": "black",  "path": f"data"},
    {"label": "Low-Fidelity", "ebsd_id": "ebsd_4", "colour": "tab:red",  "path": f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)/250318014435_i21_simulate"},
    # {"label": "High-Fidelity", "ebsd_id": "ebsd_2", "colour": "tab:blue", "path": f"{MOOSE_PATH}/2025-03-15 (617_s3_vh_x_hr)"},

    # LH2 Model
    # {"label": "Run 1", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-03-25 (lh2_x_sm8_i19)/250323214745_i7_simulate"},
    # {"label": "Run 2", "ebsd_id": "ebsd_4", "colour": "tab:red",   "path": f"{ASMBO_PATH}/2025-03-28 (lh2_x_sm8_i29)/250327093649_i16_simulate"}, # OPT
    # {"label": "Run 3", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-03-31 (lh2_x_sm8_i31)/250330063457_i18_simulate"},
    # {"label": "Run 4", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-04-02 (lh2_x_sm8_i23)/250401051359_i10_simulate"},
    # {"label": "Run 5", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-03-31 (lh2_x_sm8_i31)/250330213453_i29_simulate"},
    # {"label": "Low-Fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:red",  "path": f"{ASMBO_PATH}/2025-03-28 (lh2_x_sm8_i29)/250327093649_i16_simulate"},
    # {"label": "High-Fidelity", "ebsd_id": "ebsd_2", "colour": "tab:blue", "path": f"{MOOSE_PATH}/2025-04-05 (617_s3_lh2_di_x_hr)"},

    # LH6 Model
    # {"label": "Run 1", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-04-09 (lh6_x_sm8_i44)/250407052902_i6_simulate"},
    # {"label": "Run 2", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-04-14 (lh6_x_sm8_i32)/250413031321_i23_simulate"},
    # {"label": "Run 3", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-04-18 (lh6_x_sm8_i27)/250418123844_i27_simulate"},
    # {"label": "Run 4", "ebsd_id": "ebsd_4", "colour": "tab:red",   "path": f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250422034348_i36_simulate"},
    # {"label": "Run 5", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250420224600_i20_simulate"},
    # {"label": "Low-Fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:red",  "path": f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250422034348_i36_simulate"},
    # {"label": "High-Fidelity", "ebsd_id": "ebsd_2", "colour": "tab:blue", "path": f"{MOOSE_PATH}/2025-04-28 (617_s3_lh6_di_x_hr)"},

    # # All Models (Low-Fidelity)
    # {"label": "VH",  "ebsd_id": "ebsd_4", "colour": "tab:cyan",   "path": f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)/250318014435_i21_simulate"},
    # {"label": "LH2", "ebsd_id": "ebsd_4", "colour": "tab:orange", "path": f"{ASMBO_PATH}/2025-03-28 (lh2_x_sm8_i29)/250327093649_i16_simulate"},
    # {"label": "LH6", "ebsd_id": "ebsd_4", "colour": "tab:purple", "path": f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250422034348_i36_simulate"},

    # VH Model, cross-validation 1
    # {"label": "Run 1", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-06-03 (vh_x_sm8_i97_val)/250601025336_i19_simulate"},
    # {"label": "Run 2", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i13_val)/250609093847_i7_simulate"},
    # {"label": "Run 3", "ebsd_id": "ebsd_4", "colour": "tab:red",   "path": f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310145710_i22_simulate"},
    # {"label": "Run 4", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i52_val)/250608004520_i6_simulate"},
    # {"label": "Run 5", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i52_val)/250609222901_i8_simulate"}, #
    # # {"label": "Low-Fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i52_val)/250609222901_i8_simulate"},
    # # {"label": "High-Fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/2025-06-12 (617_s3_vh_di_x_hr_val)"},

    # VH Model, cross-validation 2
    # {"label": "Run 1", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i14_cv2)/250708043252_i4_simulate"},
    # {"label": "Run 2", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i26_cv2)/250707055737_i4_simulate"},
    # {"label": "Run 3", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i29_cv2)/250707093811_i6_simulate"},
    # {"label": "Run 4", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i30_cv2)/250708071533_i13_simulate"}, #
    # {"label": "Run 5", "ebsd_id": "ebsd_4", "colour": "tab:red",   "path": f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i32_cv2)/250708043852_i9_simulate"},
    # {"label": "Low-Fidelity",  "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i30_cv2)/250708071533_i13_simulate"},
    # {"label": "High-Fidelity", "ebsd_id": "ebsd_2", "colour": "tab:red",   "path": f"{MOOSE_PATH}/2025-07-10 (617_s3_vh_di_x_hr_val2)"},
    
    # VH Model, cross-validation 2 (SS)
    # {"label": "Run 1", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i52_val)/250609222901_i8_simulate"},
    # {"label": "Run 2", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310145710_i22_simulate"},
    # {"label": "Run 3", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310161708_i25_simulate"},
    # {"label": "Run 4", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)/250318014435_i21_simulate"},
    # {"label": "Run 5", "ebsd_id": "ebsd_4", "colour": "tab:red",   "path": f"{ASMBO_PATH}/2025-03-25 (vh_x_sm8_i31)/250325072901_i16_simulate"},
    # {"label": "Low-Fidelity", "ebsd_id": "ebsd_4", "colour": "tab:green", "path": f"{ASMBO_PATH}/2025-03-25 (vh_x_sm8_i31)/250325072901_i16_simulate"},
    # {"label": "High-Fidelity", "ebsd_id": "ebsd_4", "colour": "tab:red", "path": f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310161708_i25_simulate"},
]
for si in SIM_INFO_LIST:
    si["data"] = csv_to_dict(f"{si['path']}/summary.csv")

# Grain IDs
GRAIN_IDS = [
    [14, 72, 95, 101, 207, 240, 262, 287],  # Calibration
    [39, 50, 138, 164, 185, 223, 243, 238], # Validation
]

# Plotting parameters
OPT_ALPHA = 1.0
OTH_ALPHA = 0.4
# SPACING, FORCE_OPT, KNEE_POINT = -5.25, [0,1,2,3,4], None # Runs
SPACING, FORCE_OPT, KNEE_POINT = -2.25, [0,1], -1 # Fidelity
# SPACING, FORCE_OPT, KNEE_POINT = -6.25, [0,1,2], -1 # Models

# Script parameters
SHOW_LEGEND   = True
SHOW_GRAIN_ID = False
SHOW_ERROR    = True
SHOW_OTH      = True # whether to show optimal keys in the legend
STRAIN_FIELD  = "average_strain"
STRESS_FIELD  = "average_stress"
RES_DATA_MAP  = "data/res_grain_map.csv"

# Main function
def main():

    # Prepare experimental data
    exp_dict = csv_to_dict(EXP_PATH)
    eval_strains = np.linspace(0, exp_dict["strain_intervals"][-1], 32)

    # Calculate simulation errors
    ge_list = get_geodesic_errors(SIM_INFO_LIST, exp_dict, eval_strains, GRAIN_IDS[0]) # calibration
    se_list = get_stress_errors([si["data"] for si in SIM_INFO_LIST], exp_dict, eval_strains)
    # se_list[0] *= 1.1

    # Calculate knee point
    if KNEE_POINT != None:
        kp_index = KNEE_POINT
    else:
        kp_index = find_knee_point([ge_list, ge_list, se_list])
        print(f"Knee Point at: Run {kp_index+1}")
    
    # Determine alpha list
    alpha_list = [OTH_ALPHA]*kp_index + [OPT_ALPHA] + [OTH_ALPHA]*(len(SIM_INFO_LIST)-kp_index-1)
    if FORCE_OPT != None:
        for i in FORCE_OPT:
            alpha_list[i] = OPT_ALPHA

    # Initialise stress-strain plot
    plotter = Plotter("strain", "stress", "mm/mm", "MPa")
    plotter.prep_plot(size=14)
    plotter.set_limits((0,0.5), (0,1800))

    # Plot stress-strain data
    plt.scatter(exp_dict["strain"], exp_dict["stress"], color=EXP_COLOUR, s=8**2)
    for j, (si, alpha) in enumerate(zip(SIM_INFO_LIST, alpha_list)):
        zorder = 4 if j == kp_index else 3
        plt.plot(si["data"][STRAIN_FIELD], si["data"][STRESS_FIELD], color=si["colour"], alpha=alpha, linewidth=3, zorder=zorder)

    # Add stress errors to legend
    if SHOW_LEGEND:
        se_list = [f"{round(se*100, 1)}" for se in se_list]
        add_supp_legend(se_list, alpha_list, SPACING, "%")

    # Format and save
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    save_plot("results/plot_opts_ss.png")

    # Plot reorientation trajectories
    for i, grain_ids in enumerate(GRAIN_IDS):

        # Initialise IPF
        ipf = IPF(get_lattice("fcc"))
        direction = [1,0,0]
        get_trajectories = lambda dict, g_ids : [transpose([dict[f"g{g_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]) for g_id in g_ids]

        # Plot experimental reorientation trajectories
        exp_trajectories = get_trajectories(exp_dict, grain_ids)
        if 14 in grain_ids:
            exp_trajectories[grain_ids.index(14)] = exp_trajectories[grain_ids.index(14)][:-1] 
        if 287 in grain_ids:
            exp_trajectories[grain_ids.index(287)] = exp_trajectories[grain_ids.index(287)][:-3] 
        ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": EXP_COLOUR, "linewidth": 3})
        ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": EXP_COLOUR, "head_width": 0.01, "head_length": 0.015})
        ipf.plot_ipf_trajectory([[et[0]] for et in exp_trajectories], direction, "scatter", {"color": EXP_COLOUR, "s": 8**2})

        # Iterate through simulations
        for j, (si, alpha) in enumerate(zip(SIM_INFO_LIST, alpha_list)):

            # Get simulated reorientation trajectories
            sim_grain_ids = [get_sim_grain_id(grain_id, si["ebsd_id"]) for grain_id in grain_ids]
            sim_trajectories = get_trajectories(si["data"], sim_grain_ids)
            if 44 in sim_grain_ids:
                sim_trajectories[sim_grain_ids.index(44)] = sim_trajectories[sim_grain_ids.index(44)][:-9] 
            
            # Determine zorder
            zorder = 4 if j == kp_index else 3

            # Plot simulated reorientation trajectories
            sim_colour = si["colour"]
            ipf.plot_ipf_trajectory(sim_trajectories, direction, "plot", {"color": sim_colour, "linewidth": 2, "zorder": zorder, "alpha": alpha})
            ipf.plot_ipf_trajectory(sim_trajectories, direction, "arrow", {"color": sim_colour, "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": zorder, "alpha": alpha})
            ipf.plot_ipf_trajectory([[st[0]] for st in sim_trajectories], direction, "scatter", {"color": sim_colour, "s": 6**2, "zorder": zorder, "alpha": alpha})

            # Plot grain IDs
            if SHOW_GRAIN_ID:
                for exp_trajectory, grain_id in zip(exp_trajectories, grain_ids):
                    ipf.plot_ipf_trajectory([[exp_trajectory[0]]], direction, "text", {"color": "blue", "fontsize": 8, "s": grain_id, "zorder": 3})

        # Add geodesic errors to legend
        if SHOW_LEGEND:
            ge_list = get_geodesic_errors(SIM_INFO_LIST, exp_dict, eval_strains, grain_ids)
            ge_list = [f"{round(ge*180/np.pi, 1)}" for ge in ge_list]
            add_supp_legend(ge_list, alpha_list, SPACING, "°")
        
        # Save IPF plot
        save_plot(f"results/plot_opts_rt_{i+1}.png")

def add_supp_legend(error_list:list, alpha_list:list, spacing:float=-5.5, units:str="") -> None:
    """
    Adds supplementary information to a legend that is aligned
    horizontally to the main keys; note that this immediately
    adds the legend to the current axis

    Parameters:
    * `error_list`: Error values to add as supplementary information
    * `alpha_list`: List of alpha values for each simulation
    * `spacing`:    Spacing between main keys and information
    * `units`:      Units for the errors
    """

    # Define main keys of the legend
    handles = [plt.scatter([], [], color=EXP_COLOUR, label="Experiment", s=8**2)]
    if SHOW_OTH:
        handles += [plt.plot([], [], color=si["colour"], label=si['label'], alpha=alpha, linewidth=3)[0] for si, alpha in zip(SIM_INFO_LIST, alpha_list)]
    else:
        handles += [plt.plot([], [], color=si["colour"], label=si['label'], alpha=alpha, linewidth=3)[0] for si, alpha in zip(SIM_INFO_LIST, alpha_list)
                    if alpha == OPT_ALPHA]

    # Define supplementary information
    if SHOW_ERROR:
        formatted_error_list = [f"{error}{units}" for error in error_list]
        if not SHOW_OTH:
            formatted_error_list = [fe for fe, alpha in zip(formatted_error_list, alpha_list) if alpha == OPT_ALPHA]
        se_label_list = [" "] + [f"({fe})" for fe in formatted_error_list]
        handles += [plt.scatter([], [], color="white", label=se_label, marker="o", s=0) for se_label in se_label_list]

    # Add legend
    num_cols = 2 if SHOW_ERROR else 1
    legend = plt.legend(handles=handles, ncol=num_cols, columnspacing=spacing, framealpha=1, edgecolor="black",
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

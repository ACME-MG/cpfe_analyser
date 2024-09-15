"""
 Title:         Analyse
 Description:   Analyses the simulation results
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += ["/home/janzen/code/deer_sim"]
import numpy as np
import matplotlib.pyplot as plt
from deer_sim.analyse.plotter import Plotter, save_plot, define_legend
from deer_sim.analyse.pole_figure import PF, IPF, get_lattice, get_colour_map
from deer_sim.maths.familiser import get_grain_family
from deer_sim.helper.general import transpose, remove_consecutive_duplicates
from deer_sim.helper.io import csv_to_dict, safe_mkdir

# Constants
EXP_PATH = "../data/617_s3/617_s3_exp.csv"
MAP_PATH = "../data/617_s3/grain_map.csv"
# SIM_PATH = "./data/summary_617_s3.csv"
SIM_PATH = "./data/summary.csv"
EVP_PATH = "./data/evp_data.csv"

def get_grain_ids(exp_path:str, mesh_path:str) -> dict:
    """
    Gets the mappable grain IDs
    
    Parameters:
    * `exp_path`:  Path to the experimental data
    * `mesh_path`: Path to the mesh CSV map

    Returns a dictionary mapping the experimental grain IDs to the mesh grain IDs
    """

    # Read files
    exp_dict  = csv_to_dict(exp_path)  # exp
    mesh_dict = csv_to_dict(mesh_path) # exp : mesh
    exp_grain_ids = [int(key.replace("_phi_1","").replace("g","")) for key in exp_dict.keys() if "_phi_1" in key]
    
    # Map experimental grain IDs to mesh grain IDs
    exp_to_mesh = {}
    for exp_grain_id in exp_grain_ids:
        if exp_grain_id in mesh_dict["ebsd_id"]:
            ebsd_index = mesh_dict["ebsd_id"].index(exp_grain_id)
            exp_to_mesh[exp_grain_id] = int(mesh_dict["mesh_id"][ebsd_index])

    # Return mapping
    return exp_to_mesh

def get_trajectories(data_dict:dict, include:list=None) -> list:
    """
    Gets the reorientation trajectories

    Parameters:
    * `data_dict`: The data dictionary
    * `include`:   List of grain IDs to include;
                   uses all grain IDs if parameter unspecified
    
    Return trajectories as a list of lists of euler angles
    """

    # Read data and specify grain IDs if defined
    grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in data_dict.keys() if "_phi_1" in key]
    if include != None:
        grain_ids = [grain_id for grain_id in grain_ids if grain_id in include]

    # Get trajectories
    trajectories = []
    for grain_id in grain_ids:
        trajectory = [data_dict[f"g{grain_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]
        trajectory = transpose(trajectory)
        trajectory = remove_consecutive_duplicates(trajectory)
        trajectories.append(trajectory)

    # Return trajectories
    return trajectories

def save_plot_results(file_path:str, dir_path:str="results") -> None:
    """
    Saves the plots to a specific folder

    Parameters:
    * `file_path`: The path to the file
    * `dir_path`:  The path to the directory
    """
    safe_mkdir(dir_path)
    save_plot(f"./{dir_path}/{file_path}")

# Define specific grain IDs
grain_id_dict = get_grain_ids(EXP_PATH, MAP_PATH)
exp_grain_ids  = [159, 166, 167, 173, 180]
# exp_grain_ids  = [16, 21, 37, 46, 76, 82, 87, 99, 101, 110, 137, 141, 147, 152, 154, 159, 166, 167, 173, 180] # GOOD
# exp_grain_ids  = [23, 27, 36, 38, 40, 49, 56, 64, 66, 97, 108, 109, 112, 114, 120, 128, 130, 139, 148, 176, 178] # OKAY
sim_grain_ids  = [grain_id_dict[exp_grain_id] for exp_grain_id in exp_grain_ids]
# [print(f" exp: {exp}\t sim: {sim}") for exp, sim in zip(exp_grain_ids, sim_grain_ids)]

# Get experimental data
exp_dict = csv_to_dict(EXP_PATH)
exp_ss   = {"strain": exp_dict["strain"], "stress": exp_dict["stress"]}
exp_trajectories = get_trajectories(exp_dict, exp_grain_ids)
exp_start_orientations = [trajectory[0] for trajectory in get_trajectories(exp_dict)]
exp_final_orientations = [trajectory[-1] for trajectory in get_trajectories(exp_dict)]

# Get simulated data
sim_dict     = csv_to_dict(SIM_PATH)
sim_grain_ss = {"strain": sim_dict["average_grain_strain"], "stress": sim_dict["average_grain_stress"]}
sim_grip_ss  = {"strain": sim_dict["average_grip_strain"], "stress": sim_dict["average_grip_stress"]}
sim_trajectories       = get_trajectories(sim_dict, sim_grain_ids)
sim_start_orientations = [trajectory[0] for trajectory in get_trajectories(sim_dict)]
sim_final_orientations = [trajectory[-1] for trajectory in get_trajectories(sim_dict)]

# Gets EVP data
evp_dict = csv_to_dict(EVP_PATH)
evp_strain = [strain for strain in evp_dict["cal_tensile_strain"] if strain < 0.45]
evp_stress = evp_dict["cal_tensile_stress"][:len(evp_strain)]
evp_ss = {"strain": evp_strain, "stress": evp_stress}

# Initialise IPF and PF plotters
ipf = IPF(get_lattice("fcc"))
pf = PF(get_lattice("fcc"))
direction = [1,0,0]

# Plot stress-strain curves
plotter_ss = Plotter("strain", "stress", "mm/mm", "MPa")
plotter_ss.prep_plot()
plotter_ss.scat_plot(exp_ss,       colour="darkgray", name="Experimental")
plotter_ss.line_plot(evp_ss,       colour="green",    name="VP (Taylor)")
plotter_ss.line_plot(sim_grain_ss, colour="red",      name="CP (FEM)")
plotter_ss.line_plot(sim_grip_ss,  colour="blue",     name="VP (FEM)")
plotter_ss.set_legend()
save_plot_results("plot_ss.png")

# Plot PF
for plane in [[1,0,0], [1,1,0], [1,1,1]]:
    plane_str = "".join([str(p) for p in plane])
    pf.plot_pf(exp_final_orientations, plane)
    save_plot_results(f"plot_pf_exp_{plane_str}.png")
    pf.plot_pf(sim_final_orientations, plane)
    save_plot_results(f"plot_pf_sim_{plane_str}.png")

# Plot final stress distribution on IPF
sim_final_stresses = [sim_dict[key][-1] for key in sim_dict.keys() if key.startswith("g") and "_stress" in key]
ipf.plot_ipf(sim_final_orientations, direction, sim_final_stresses)
save_plot_results("plot_ipf_stress.png")
get_colour_map(sim_final_stresses, "vertical")
save_plot_results("plot_ipf_stress_cm.png")

# Plot initial orientations on IPF
ipf.plot_ipf_trajectory([[eso] for eso in exp_start_orientations], direction, "scatter", {"color": "darkgray", "s": 8**2})
ipf.plot_ipf_trajectory([[sso] for sso in sim_start_orientations], direction, "scatter", {"color": "green", "s": 6**2, "zorder": 3})
define_legend(["darkgray", "green"], ["Experimental", "CP (FEM)"], ["scatter", "scatter"])
save_plot_results("plot_ipf_initial.png")

# Plot final orientations on IPF
ipf.plot_ipf_trajectory([[efo] for efo in exp_final_orientations], direction, "scatter", {"color": "darkgray", "s": 8**2})
ipf.plot_ipf_trajectory([[sfo] for sfo in sim_final_orientations], direction, "scatter", {"color": "green", "s": 6**2, "zorder": 3})
define_legend(["darkgray", "green"], ["Experimental", "CP (FEM)"], ["scatter", "scatter"])
save_plot_results("plot_ipf_final.png")

# Plot trajectories
ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": "darkgray", "linewidth": 2})
ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": "darkgray", "head_width": 0.01, "head_length": 0.015})
ipf.plot_ipf_trajectory([[et[0]] for et in exp_trajectories], direction, "scatter", {"color": "darkgray", "s": 8**2})
for exp_trajectory, grain_id in zip(exp_trajectories, exp_grain_ids):
    ipf.plot_ipf_trajectory([[exp_trajectory[0]]], direction, "text", {"color": "black", "fontsize": 8, "s": grain_id})
ipf.plot_ipf_trajectory(sim_trajectories, direction, "plot", {"color": "green", "linewidth": 1, "zorder": 3})
ipf.plot_ipf_trajectory(sim_trajectories, direction, "arrow", {"color": "green", "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
ipf.plot_ipf_trajectory([[st[0]] for st in sim_trajectories], direction, "scatter", {"color": "green", "s": 6**2, "zorder": 3})
define_legend(["darkgray", "green"], ["Experimental", "CP (FEM)"], ["line", "line"])
save_plot_results("plot_ipf_trajectories.png")

# Initialise elastic strain / stress plotting
crystal_directions = [[2,2,0], [1,1,1], [3,1,1], [2,0,0], [3,3,1], [5,3,1]]
colour_list = ["green", "black", "blue", "red", "purple", "orange"]
plotter_es = Plotter("Elastic Strain", "Applied Stress", "µmm/mm", "MPa")
plotter_es.prep_plot()
sim_elastics = [sim_dict[key] for key in sim_dict.keys() if key.startswith("g") and "_elastic" in key]
sim_volumes = [sim_dict[key] for key in sim_dict.keys() if key.startswith("g") and "_volume" in key]

# Plot elastic strains and stresses
for crystal_direction, colour in zip(crystal_directions, colour_list):
    family_indices = get_grain_family(sim_start_orientations, crystal_direction, direction, 10)
    family_elastics = transpose([sim_elastics[i] for i in family_indices])
    family_volumes = transpose([sim_volumes[i] for i in family_indices])
    average_elastics = [np.average(family_elastic, weights=family_volume) if sum(family_volume) > 0 else np.average(family_elastic)
                        for family_elastic, family_volume in zip(family_elastics, family_volumes)]
    if len(average_elastics) == 0:
        continue
    average_elastics = [ae*1e6 for ae in average_elastics]
    average_dict = {"Elastic Strain": average_elastics, "Applied Stress": sim_dict["average_grain_stress"]}
    crystal_str = "{" + "".join([str(cd) for cd in crystal_direction]) + "}"
    plotter_es.scat_plot(average_dict, colour=colour, name=crystal_str)
    plotter_es.line_plot(average_dict, colour=colour)
plotter_es.set_legend()
save_plot_results("plot_es.png")

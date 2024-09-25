"""
 Title:         Elastic
 Description:   Plots the elastic strain
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += ["/home/janzen/code/moose_sim"]
import numpy as np
from moose_sim.helper.general import transpose
from moose_sim.helper.io import csv_to_dict
from moose_sim.analyse.plotter import Plotter, save_plot
from moose_sim.maths.familiser import get_grain_family
from moose_sim.analyse.pole_figure import IPF, get_lattice, get_colour_map

# Read data
SUMMARY_PATH = "data/summary_ae.csv"
summary_dict = csv_to_dict(SUMMARY_PATH)

# Get orientations (euler-bunge, passive, rads)
grain_ids = [int(key.replace("_phi_1","").replace("g","")) for key in summary_dict.keys() if "_phi_1" in key]
orientation_keys = [[f"g{grain_id}_{phi}" for phi in ["phi_1", "Phi", "phi_2"]] for grain_id in grain_ids]
start_orientations = [[summary_dict[key][0] for key in keys] for keys in orientation_keys]
final_orientations = [[summary_dict[key][-1] for key in keys] for keys in orientation_keys]

# Get elastic strain
elastics = [summary_dict[key] for key in summary_dict.keys() if key.startswith("g") and "_elastic" in key]
volumes = [summary_dict[key] for key in summary_dict.keys() if key.startswith("g") and "_volume" in key]

# Initialise elastic strain / stress plotting
sample_direction = [1,0,0]
crystal_directions = [[2,2,0], [1,1,1], [3,1,1], [2,0,0]]
colour_list = ["green", "black", "blue", "red"]
plotter_es = Plotter("Elastic Strain", "Applied Stress", "mm/mm", "MPa")
plotter_es.prep_plot()

# Plot elastic strains and stresses
for crystal_direction, colour in zip(crystal_directions, colour_list):
    family_indices = get_grain_family(start_orientations, crystal_direction, sample_direction, 10)
    family_elastics = transpose([elastics[i] for i in family_indices])
    family_volumes = transpose([volumes[i] for i in family_indices])
    average_elastics = [np.average(family_elastic, weights=family_volume) if sum(family_volume) > 0 else np.average(family_elastic)
                        for family_elastic, family_volume in zip(family_elastics, family_volumes)]
    average_dict = {"Elastic Strain": average_elastics, "Applied Stress": summary_dict["average_grain_stress"]}
    crystal_str = "{" + "".join([str(cd) for cd in crystal_direction]) + "}"
    plotter_es.scat_plot(average_dict, colour=colour, name=crystal_str)
    plotter_es.line_plot(average_dict, colour=colour)
plotter_es.set_legend()
save_plot("results/plot_es.png")
    
# Plot stress distribution
ipf = IPF(get_lattice("fcc"))
final_stresses = [summary_dict[key][-1] for key in summary_dict.keys() if key.startswith("g") and "_stress" in key]
ipf.plot_ipf(final_orientations, sample_direction, final_stresses)
save_plot("results/plot_ipf_stress.png")
get_colour_map(final_stresses, "vertical")
save_plot("results/plot_ipf_stress_cm.png")

# Plot stress-strain
grain_ss = {"strain": summary_dict["average_grain_strain"], "stress": summary_dict["average_grain_stress"]}
grip_ss  = {"strain": summary_dict["average_grip_strain"], "stress": summary_dict["average_grip_stress"]}
plotter_ss = Plotter("strain", "stress", "mm/mm", "MPa")
plotter_ss.prep_plot()
plotter_ss.line_plot(grain_ss, colour="red",  name="CP (FEM)")
plotter_ss.line_plot(grip_ss,  colour="blue", name="VP (FEM)")
plotter_ss.set_legend()
save_plot("results/plot_ss.png")
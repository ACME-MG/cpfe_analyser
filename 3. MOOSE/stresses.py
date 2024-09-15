"""
 Title:         Elastic
 Description:   Plots the elastic strain
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += ["/home/janzen/code/deer_sim"]
import cv2, numpy as np
import matplotlib.pyplot as plt
from deer_sim.helper.general import flatten
from deer_sim.helper.io import csv_to_dict
from deer_sim.analyse.plotter import save_plot
from deer_sim.analyse.pole_figure import IPF, get_lattice, get_colour_map

# Read data
SUMMARY_PATH = "data/sim_data.csv"
KEYWORD      = "_stress"
summary_dict = csv_to_dict(SUMMARY_PATH)

# Extract initialisation information
grain_ids = [int(key.replace("_phi_1","").replace("g","")) for key in summary_dict.keys() if "_phi_1" in key]
orientation_keys = [[f"g{grain_id}_{phi}" for phi in ["phi_1", "Phi", "phi_2"]] for grain_id in grain_ids]
num_states = len(summary_dict[orientation_keys[0][0]])

# Get orientations (euler-bunge, passive, rads) and grain stress for each state
orientation_history = [[[summary_dict[key][i] for key in keys] for keys in orientation_keys]
                       for i in range(num_states)]
stress_history = [[summary_dict[key][i] for key in summary_dict.keys()
                   if key.startswith("g") and KEYWORD in key] for i in range(num_states)]
volume_history = [[summary_dict[key][i] for key in summary_dict.keys()
                   if key.startswith("g") and "volume" in key] for i in range(num_states)]
all_stresses = flatten(stress_history)
all_volumes = flatten(volume_history)

# Initialise video writer
frame_rate = 3
frame_size = (640, 480)
fourcc = cv2.VideoWriter_fourcc(*"mp4v") # codec for mp4
video_writer = cv2.VideoWriter("stresses.mp4", fourcc, frame_rate, frame_size)

# Record stress distribution
ipf = IPF(
    lattice       = get_lattice("fcc"),
    colour_limits = (min(all_stresses), max(all_stresses)),
    size_limits   = (min(all_volumes), max(all_volumes)),
)
sample_direction = [1,0,0]
for i, (orientations, stresses, volumes) in enumerate(zip(orientation_history, stress_history, volume_history)):
    
    # Ignore initial (zero stress)
    if i == 0:
        continue
    print(f"Adding frame {i}")
    
    # Otherwise, draw the IPF
    ipf.plot_ipf(orientations, sample_direction, stresses, volumes)
    figure = plt.gcf()
    figure.canvas.draw()
    
    # Convert drawing to video frame
    image = np.frombuffer(figure.canvas.tostring_rgb(), dtype=np.uint8)
    image = image.reshape(figure.canvas.get_width_height()[::-1] + (3,))
    image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)
    video_writer.write(image)
    plt.close(figure)
    
# Release writer
video_writer.release()

get_colour_map(all_stresses, orientation="vertical")
save_plot("results/cm.png")

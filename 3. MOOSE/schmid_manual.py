"""
 Title:         Elastic
 Description:   Plots the elastic strain
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += ["/home/janzen/code/moose_sim"]
import cv2, math, numpy as np
import matplotlib.pyplot as plt
from moose_sim.maths.neml import reorient_euler
from moose_sim.helper.general import flatten
from moose_sim.helper.io import csv_to_dict
from moose_sim.analyse.plotter import save_plot
from moose_sim.analyse.pole_figure import IPF, get_lattice, get_colour_map
from neml.math import rotations, tensors
from neml.cp import crystallography

# Constants
SUMMARY_PATH = "data/sim_data.csv"
KEYWORD      = "_stress"

def euler_to_matrix(euler:list) -> list:
    """
    Determines the orientation matrix of a set of euler-bunge angles (rads);
    from Ref. [1]
    
    Parameters:
    * `euler`: The euler angle in euler-bunge form

    Returns the 3x3 orientation matrix as a list of lists
    """
    om_11 = math.cos(euler[0])*math.cos(euler[2]) - math.sin(euler[0])*math.sin(euler[2])*math.cos(euler[1])
    om_12 = math.sin(euler[0])*math.cos(euler[2]) + math.cos(euler[0])*math.sin(euler[2])*math.cos(euler[1])
    om_13 = math.sin(euler[2])*math.sin(euler[1])
    om_21 = -math.cos(euler[0])*math.sin(euler[2]) - math.sin(euler[0])*math.cos(euler[2])*math.cos(euler[1])
    om_22 = -math.sin(euler[0])*math.sin(euler[2]) + math.cos(euler[0])*math.cos(euler[2])*math.cos(euler[1])
    om_23 = math.cos(euler[2])*math.sin(euler[1])
    om_31 = math.sin(euler[0])*math.sin(euler[1])
    om_32 = -math.cos(euler[0])*math.sin(euler[1])
    om_33 = math.cos(euler[1])
    om = [[om_11, om_12, om_13],
          [om_21, om_22, om_23],
          [om_31, om_32, om_33]]
    return om

def get_lattice(structure:str="fcc"):
    """
    Gets the lattice object

    Parameters:
    * `structure`: The crystal structure

    Returns the lattice object
    """
    lattice = crystallography.CubicLattice(1.0)
    if structure == "fcc":
        lattice.add_slip_system([1,1,0], [1,1,1])
    elif structure == "bcc":
        lattice.add_slip_system([1,1,1], [1,1,0])
        lattice.add_slip_system([1,1,1], [1,2,3])
        lattice.add_slip_system([1,1,1], [1,1,2])
    else:
        raise ValueError(f"Crystal structure '{structure}' unsupported!")
    return lattice

def get_schmid(euler:list) -> float:
    """
    Gets the Schmid factor
    
    Parameters:
    * `euler`: The euler-bunge angle (rads)
    
    Returns the Schmid factor
    """
    
    # Iniitalise
    get_term = lambda u_1, v_1, w_1, u_2, v_2, w_2 : (u_1*u_2+v_1*v_2+w_1*w_2)/math.sqrt((u_1**2+v_1**2+w_1**2)*(u_2**2+v_2**2+w_2**2))
    matrix = euler_to_matrix(euler)
    vector = np.dot(np.array(matrix), np.array(direction))
    
    # Calculate schmid factors
    schmid_factors = []
    for fcc_plane in fcc_planes:
        term_1 = get_term(*fcc_plane, *vector)
        for fcc_direction in fcc_directions:
            term_2 = get_term(*vector, *fcc_direction)
            schmid_factor = term_1*term_2
            schmid_factors.append(schmid_factor)
    
    # Get averaged schmid factors and return
    average_schmid = np.average(schmid_factors)
    return average_schmid
    
# Extract initialisation information
summary_dict = csv_to_dict(SUMMARY_PATH)
grain_ids = [int(key.replace("_phi_1","").replace("g","")) for key in summary_dict.keys() if "_phi_1" in key]
orientation_keys = [[f"g{grain_id}_{phi}" for phi in ["phi_1", "Phi", "phi_2"]] for grain_id in grain_ids]
num_states = len(summary_dict[orientation_keys[0][0]])

# Get orientations (euler-bunge, passive, rads) and volumes at each state
orientation_history = [[[summary_dict[key][i] for key in keys] for keys in orientation_keys]
                       for i in range(num_states)]
volume_history = [[summary_dict[key][i] for key in summary_dict.keys()
                   if key.startswith("g") and "volume" in key] for i in range(num_states)]
all_volumes = flatten(volume_history)

# Calculate schmid factor history
lattice = get_lattice("fcc")
direction = [1,0,0]
schmid_history = []
# fcc_directions = [[1,-1,0], [0,1,-1], [1,0,-1], [1,1,0], [0,1,1], [1,0,1]]
fcc_planes = [[1,1,1], [1,-1,1], [-1,1,1], [-1,-1,1]]
fcc_directions = [[1,-1,0], [0,1,-1], [1,0,-1]]

# Iterate through orientations
for orientation_list in orientation_history:
    schmid_list = []
    for orientation in orientation_list:
        euler = reorient_euler(orientation)
        schmid = get_schmid(euler)
        schmid_list.append(schmid)
    schmid_history.append(schmid_list)
all_schmids = flatten(schmid_history)

# Initialise video writer
frame_rate = 3
frame_size = (640, 480)
fourcc = cv2.VideoWriter_fourcc(*"mp4v") # codec for mp4
video_writer = cv2.VideoWriter("schmid.mp4", fourcc, frame_rate, frame_size)

# Record stress distribution
ipf = IPF(
    lattice = get_lattice("fcc"),
    # colour_limits = (min(all_schmids), max(all_schmids)),
    # size_limits = (min(all_volumes), max(all_volumes)),
)
sample_direction = [1,0,0]
for i, (orientations, schmids, volumes) in enumerate(zip(orientation_history, schmid_history, volume_history)):
    
    # Ignore initial (zero stress)
    if i == 0:
        continue
    print(f"Adding frame {i}")
    
    # Otherwise, draw the IPF
    # ipf.plot_ipf(orientations, sample_direction, schmids, volumes)
    ipf.plot_ipf(orientations, sample_direction, schmids)
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

get_colour_map(all_schmids, orientation="vertical")
# save_plot("results/cm.png")

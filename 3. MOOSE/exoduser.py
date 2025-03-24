"""
 Title:         Exoduser
 Description:   Plots the simulated deformation using an exodus file
 Author:        Janzen Choi

"""

# Libraries
import itertools
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import pyvista as pv
import time
from copy import deepcopy
import sys; sys.path += [".."]
from __common__.general import integer_to_ordinal
from __common__.io import csv_to_dict
from __common__.ipf_cubic import euler_to_rgb
from __common__.orientation import deg_to_rad, rad_to_deg
from __common__.plotter import save_plot
from __common__.video import Video

# Resolution information
MESH_INDEX = 1
MESH_INFO = [

    # Low-fidelity mesh
    {
        "exodus":    "data/617_s3/40um/mesh.e",
        "path":      "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/data/2024-06-26 (ansto_617_s3)/prior_with_stage/res20gs5/ebsdExportColumnsTableReduced_FillRegion.csv",
        "step":      20,
        "thickness": 3,
        "dec_res":   2,
    },

    # High-fidelity mesh
    {
        "exodus":    "data/617_s3/10um/mesh.e",
        "path":      "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/data/2024-06-26 (ansto_617_s3)/prior_with_stage/res10gs10/ebsdExportColumnsTableReduced_FillRegion.csv",
        "step":      10,
        "thickness": 10,
        "dec_res":   1,
    },
][MESH_INDEX]

# Simulation paths
ASMBO_PATH    = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo"
MOOSE_PATH    = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim"
# SIM_PATH      = f"{MOOSE_PATH}/2025-03-19 (617_s3_vh_di_x)"
SIM_PATH      = f"{MOOSE_PATH}/2025-03-23 (617_s3_vh_di_x_hr)"
EXODUS_PREFIX = "simulation_exodus_ts"
# for file in simulation_exodus.e*; do mv "$file" "$(echo "$file" | sed -E 's/simulation_exodus\.e(-s0*([0-9]+))?/simulation_exodus_ts\2.e/')"; done

# Plotting parameters
CREATE_PLOT   = True
FIGURE_FACTOR = 40   # factor to change the size of the figure
IPF_DIRECTION = "x"  # IPF direction to colour the grains
EXODUS_FACE   = "xy" # face to be plotted (e.g., "xy")
# EXODUS_FACE   = "xz" # face to be plotted (e.g., "xy")
EXODUS_NORMAL = True # normal of the face to be plotted (True / False)

# Video parameters
CREATE_VIDEO = True
VIDEO_SIZE   = (3000, 1600)
VIDEO_LIMITS = [(0, VIDEO_SIZE[0]), (0, VIDEO_SIZE[1])]

# Script parameters (only change if using custom scripts)
EBSD_HEADERS = ["x", "y", "grainId", "Euler_phi1", "Euler_Phi", "Euler_phi2"]
EBSD_DEGREES = True
SIM_DEGREES  = False # radians
VOID_ORIENTATION = (0,0,0)
VOID_ELEMENT_ID = 100000 # large number

def main():
    """
    Main function
    """

    # Initialise
    prg = Progresser(f"{0}.")
    if CREATE_VIDEO:
        video = Video("results/mesh", 4, VIDEO_SIZE)

    # Reads the mesh
    prg.progress("Reading the initial mesh")
    init_exodus = pv.read(MESH_INFO["exodus"])[0]
    prg.end()

    # Read initial EBSD information
    prg.progress("Reading initial EBSD map")
    ebsd_step = MESH_INFO["step"]
    element_grid = read_elements(MESH_INFO["path"], ebsd_step, EBSD_HEADERS, EBSD_DEGREES)
    element_grid = decrease_resolution(element_grid, MESH_INFO["dec_res"])
    ebsd_step *= MESH_INFO["dec_res"]
    prg.end()
    
    # Check number of grains
    prg.progress("Checking grains")
    ebsd_grains = len(get_grain_ids(element_grid))
    mesh_grains = init_exodus.n_blocks
    prg.print(f"EBSD grains: {ebsd_grains}")
    prg.print(f"Mesh grains: {mesh_grains}")
    prg.end()

    # Map exodus grains to EBSD grains
    prg.progress("Mapping exodus grains to EBSD grains")
    exo_to_spn, _ = map_exo_to_spn(init_exodus, element_grid, MESH_INFO["thickness"])
    prg.print(f"Unique maps: {len(list(set(exo_to_spn.values())))}/{mesh_grains}")
    prg.end()
    
    # Map exodus elements to EBSD voxels
    prg.progress("Mapping exodus elements to EBSD voxels")
    grain_dict = get_element_info(init_exodus, element_grid, exo_to_spn, ebsd_step)
    prg.end()

    # Get elements of indexes to plot
    prg.progress("Identifying elements to plot")
    element_indexes = get_element_indexes(init_exodus, EXODUS_FACE, EXODUS_NORMAL)

    # Prepare simulation files for plotting
    prg.progress("Preparing simulation files")
    sim_dict = csv_to_dict(f"{SIM_PATH}/summary.csv")
    sim_grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in sim_dict.keys() if "_phi_1" in key]
    spn_to_exo = {v: k for k, v in exo_to_spn.items()}
    exodus_paths = [f"{SIM_PATH}/{file}" for file in os.listdir(SIM_PATH) if os.path.isfile(f"{SIM_PATH}/{file}") and EXODUS_PREFIX in file]
    exodus_paths = [f"{SIM_PATH}/{EXODUS_PREFIX}{i+1}.e" for i in range(len(exodus_paths))]
    # exodus_paths = [f"{SIM_PATH}/{EXODUS_PREFIX}{i}.e" for i in [2, 24]]
    prg.print(f"Exodus files: {len(exodus_paths)}")
    prg.print(f"Timesteps:    {len(sim_dict[list(sim_dict.keys())[0]])}")
    prg.end()

    # Iterate through timesteps
    for i, exodus_path in enumerate(exodus_paths):
        
        # Initialise
        prg = Progresser(f"{i+1}.")
        ordinal = integer_to_ordinal(i+1)

        # Apply orientations to map
        prg.progress(f"Applying orientations at timestep {i+1}")
        mappable_count = 0
        for sgi in sim_grain_ids:
            if not sgi in spn_to_exo.keys():
                continue
            euler = [sim_dict[f"g{sgi}_{phi}"][i] for phi in ["phi_1", "Phi", "phi_2"]]
            exo_id = spn_to_exo[sgi]
            for j in range(len(grain_dict[exo_id])):
                grain_dict[exo_id][j].set_orientation(*euler, SIM_DEGREES)
            mappable_count += 1
        prg.print(f"Successfully oriented {mappable_count}/{mesh_grains} grains")
        prg.end()

        # Read exodus file
        prg.progress(f"Reading the {ordinal} deformed mesh")
        prg.print(f"Located {exodus_path}")
        ts_exodus = pv.read(exodus_path)[0]
        prg.end()
        
        # Creates a plot
        if CREATE_PLOT:
            prg.progress(f"Plotting the {ordinal} deformed mesh")
            initialise_plot(ts_exodus, FIGURE_FACTOR, EXODUS_FACE)
            plot_mesh(ts_exodus, grain_dict, element_indexes, IPF_DIRECTION, EXODUS_FACE, EXODUS_NORMAL)
            plt.axis("off")
            settings = {"bbox_inches": "tight", "pad_inches": 0, "transparent": True}
            save_plot(f"results/mesh_ts{i+1}.png", settings)
            prg.end()

        # Add to video
        if CREATE_VIDEO:
            prg.progress(f"Adding the {ordinal} frame to the deformation video")
            initialise_plot(ts_exodus, FIGURE_FACTOR, EXODUS_FACE, (VIDEO_SIZE[0]/100, VIDEO_SIZE[1]/100), VIDEO_LIMITS)
            plot_mesh(ts_exodus, grain_dict, element_indexes, IPF_DIRECTION, EXODUS_FACE, EXODUS_NORMAL)
            plt.axis("off")
            figure = plt.gcf()
            video.write_to_video(figure)
            plt.cla()
            plt.clf()
            plt.close()
            prg.end()
        
# Progress updater class
class Progresser:
    def __init__(self, prefix:str=""):
        self.step = 1
        self.prefix = prefix
        self.start = time.time()
    def progress(self, message:str):
        line = "-"*20
        print(f"{line}=[ {self.prefix}{self.step} ]={line}")
        print("")
        print(f"  {message}")
        self.step += 1
    def print(self, message:str):
        print(f"  {message}")
    def end(self):
        print(f"  Finished in {round(time.time()-self.start,1)}s")
        print("")
        self.start = time.time()

def decrease_resolution(element_grid:list, factor:float) -> list:
    """
    Decreases the resolution of the voxellation
    
    Parameters:
    * `element_grid`: Grid of elements
    * `factor`:       The factor of the resolution decrease

    Returns the new element grid
    """
    new_x_size = math.ceil(len(element_grid[0]) / factor)
    new_y_size = math.ceil(len(element_grid) / factor)
    new_element_grid = get_void_element_grid(new_x_size, new_y_size)
    for row in range(new_y_size):
        for col in range(new_x_size):
            element = element_grid[math.floor(row*factor)][math.floor(col*factor)]
            new_element_grid[row][col] = deepcopy(element)
    return new_element_grid

def get_element_info(mesh, element_grid:list, exo_to_spn:dict, step_size:float) -> dict:
    """
    Gets a list of elements; the elements are ordered by the grains
    then elements within the grains
    
    Parameters:
    * `mesh`:         The mesh object
    * `element_grid`: The grid of elements
    * `exo_to_spn`:   Mapping from exodus id to spn id
    * `step_size`:    The size of each element

    Returns an ordered dict of grains containing lists of element objects
    """

    # Initialise
    position_dict = get_grain_positions(element_grid)
    get_distance = lambda a, b : math.sqrt(math.pow(a[0]-b[0],2) + math.pow(a[1]-b[1],2))
    grain_dict = {}

    # Read grains and iterate through them
    for i in range(mesh.n_blocks):

        # Get grain and grain information
        grain = mesh[i]
        exo_id = int(str(mesh.get_block_name(i)).split(" ")[-1])
        if not exo_id in exo_to_spn.keys():
            continue
        grain_id = exo_to_spn[exo_id]

        # Get centroids of the elements in the grain
        elements = grain.cell_centers().points
        centroid_list = [list(element) for element in elements]

        # Iterate through centroids
        element_list = []
        for centroid in centroid_list:

            # Get the positions of all elements in the grain
            positions = position_dict[grain_id]
            positions = [(positions["x"][i], positions["y"][i]) for i in range(len(positions["x"]))]
            positions_scaled = [(position[0]*step_size, position[1]*step_size) for position in positions]

            # Find element closest to centroid
            distances = [get_distance(position, centroid) for position in positions_scaled]
            min_index = distances.index(min(distances))
            opt_position = positions[min_index]
            element = element_grid[opt_position[1]][opt_position[0]]
            element_list.append(element)
        grain_dict[exo_id] = element_list

    # Return the dictionary of element objects
    return grain_dict

def get_grain_positions(element_grid:list) -> dict:
    """
    Gets the x and y positions for all the grains
    
    Parameters:
    * `element_grid`: A grid of elements
    
    Returns a dictionary that maps the grain IDs
    to their x and y positions
    """
    
    # Initialise dictionary to store element positions
    grain_ids = get_grain_ids(element_grid)
    position_dict = {}
    for grain_id in grain_ids:
        position_dict[grain_id] = {"x": [], "y": []}
    
    # Add points to the element dictionary
    for row in range(len(element_grid)):
        for col in range(len(element_grid[row])):
            grain_id = element_grid[row][col].get_grain_id()
            if grain_id in [VOID_ELEMENT_ID]:
                continue
            position_dict[grain_id]["x"].append(col)
            position_dict[grain_id]["y"].append(row)
    
    # Returns the position mapping
    return position_dict

def get_grain_ids(element_grid:list) -> list:
    """
    Gets the grain IDs of the element grid
    
    Parameters:
    * `element_grid`: A grid of elements
    
    Returns the list of grain IDs
    """
    grain_ids = list(set([element.get_grain_id() for element_list in element_grid for element in element_list]))
    if VOID_ELEMENT_ID in grain_ids:
        grain_ids.remove(VOID_ELEMENT_ID)
    return grain_ids

def map_exo_to_spn(mesh, element_grid:list, thickness:float) -> tuple:
    """
    Maps the grains from the Exodus file to the SPN file
    
    Parameters:
    * `mesh`:         The mesh object
    * `element_grid`: A grid of elements
    * `thickness`:    Number of elements along the z-axis
    
    Returns a mapping of the SPN to exodus IDs and the confidence of the mapping
    """

    # Reads the contents of the exodus file
    bounds = mesh.bounds
    exo_bounds = [{"min": bounds[2*i], "max": bounds[2*i+1], "range": bounds[2*i+1] - bounds[2*i]} for i in range(3)]

    # Get information from element grid
    spn_size = (len(element_grid[0]), len(element_grid), thickness)
    voxel_list = []
    for i in range(spn_size[0]):         # x
        for j in range(spn_size[1]):     # y
            for _ in range(spn_size[2]): # z
                voxel_list.append(element_grid[j][i].get_grain_id())

    # Initialise dictionaries
    exo_to_spn = {}
    confidence_list = []

    # Iterate through the exodus grains
    for i in range(mesh.n_blocks):
        
        # Get grain elements
        exo_grain = mesh[i]
        elements = exo_grain.cell_centers().points
        elements = [list(element) for element in elements]

        # Get the grain ids based on element coordinates
        id_list = []
        for element in elements:
            if math.nan in element:
                continue
            pos = [math.floor((element[j] - exo_bounds[j]["min"]) / exo_bounds[j]["range"] * spn_size[j]) for j in range(3)]
            grain_id_index = pos[0] * spn_size[1] * spn_size[2] + pos[1] * spn_size[2] + pos[2]
            grain_id = voxel_list[grain_id_index]
            id_list.append(grain_id)
        
        # Add exodus grain id
        mode = max(set(id_list), key=id_list.count)
        freq = id_list.count(mode)
        total = len(id_list)
        confidence = round(freq / total * 100, 2)

        # Update dictionaries
        exo_id = int(str(mesh.get_block_name(i)).split(" ")[-1])
        exo_to_spn[exo_id] = mode
        confidence_list.append(confidence)

    # Return
    return exo_to_spn, confidence_list

def read_elements(path:str, step_size:float, headers:list=None, degrees:bool=True) -> tuple:
    """
    Converts a CSV file into a grid of elements
    
    Parameters:
    * `path`:      The path to the CSV file
    * `step_size`: The step size of each element
    * `headers`:   The list of headers for the x coordinates, y coordinates,
                   grain IDs, phi_1, Phi, and phi_2 values
    * `degrees`:   Whether to store the orientations as degrees
    
    Returns the element grid
    """

    # Open file and read headers
    file = open(path, "r")
    csv_headers = file.readline().replace("\n", "").split(",")
    rows = file.readlines()
    
    # Get column indexes
    if headers == None:
        headers = ["x", "y", "grain_id", "phi_1", "Phi", "phi_2"]
    index_x         = csv_headers.index(headers[0])
    index_y         = csv_headers.index(headers[1])
    index_grain_id  = csv_headers.index(headers[2])
    index_avg_phi_1 = csv_headers.index(headers[3])
    index_avg_Phi   = csv_headers.index(headers[4])
    index_avg_phi_2 = csv_headers.index(headers[5])

    # Get dimensions
    x_cells, x_min = get_info([float(row.split(",")[index_x]) for row in rows], step_size)
    y_cells, y_min = get_info([float(row.split(",")[index_y]) for row in rows], step_size)
    
    # Read CSV and fill grid
    element_grid = get_void_element_grid(x_cells, y_cells)
    for row in rows:

        # Process data
        row_list = row.replace("\n", "").split(",")
        if "NaN" in row_list or "nan" in row_list:
            continue
        row_list = [float(val) for val in row_list]
        grain_id = round(row_list[index_grain_id])

        # Add to element grid
        x = round(float(row_list[index_x] - x_min) / step_size)
        y = round(float(row_list[index_y] - y_min) / step_size)
        element_grid[y][x] = Element(
            phi_1    = row_list[index_avg_phi_1],
            Phi      = row_list[index_avg_Phi],
            phi_2    = row_list[index_avg_phi_2],
            grain_id = grain_id,
            degrees  = degrees,
        )
    
    # Close file and return element grid
    file.close()
    return element_grid

def get_void_element_grid(x_cells:list, y_cells:list) -> list:
    """
    Creates a grid of void elements
    
    Parameters:
    * `x_cells`:    The number of elements on the horizontal axis
    * `y_cells`:    The number of elements on the vertical axis
    * `init_value`: The initial value of the cell in the element grid
    
    Returns a grid of void elements
    """
    element_grid = []
    for _ in range(y_cells):
        element_list = []
        for _ in range(x_cells):
            element = get_void_element()
            element_list.append(element)
        element_grid.append(element_list)
    return element_grid

def get_info(value_list:list, step_size:float) -> tuple:
    """
    Gets the range and step size from a list of values
    
    Parameters:
    * `value_list`: List of values
    * `step_size`:  The step size of each element
    
    Returns the number of values and minimum values
    """
    max_value = max(value_list)
    min_value = min(value_list)
    num_values = round((max_value - min_value) / step_size) + 1
    return num_values, min_value

# Element class
class Element:
    
    def __init__(self, phi_1:float, Phi:float, phi_2:float, grain_id:int, degrees:bool=True):
        """
        Contains information about an of a grain
        
        Parameters:
        * `phi_1`:    The phi_1 orientation of the grain
        * `Phi`:      The Phi orientation of the grain
        * `phi_2`:    The phi_2 orientation of the grain
        * `grain_id`: The ID of the grain to which the element belongs
        * `degrees`:  Whether to store the orientations as degrees
        """
        self.set_orientation(phi_1, Phi, phi_2, degrees)
        self.set_grain_id(grain_id)
    
    def set_orientation(self, phi_1:float, Phi:float, phi_2:float, degrees:bool=True) -> None:
        """
        Sets the orientations
        
        Parameters:
        * `phi_1`:   The average phi_1 orientation of the grain
        * `Phi`:     The average Phi orientation of the grain
        * `phi_2`:   The average phi_2 orientation of the grain
        * `degrees`: Whether to store the orientations as degrees
        """
        self.phi_1   = phi_1
        self.Phi     = Phi
        self.phi_2   = phi_2
        self.degrees = degrees

    def get_orientation(self, degrees:bool=True) -> tuple:
        """
        Returns the orientation of the grain
        
        Parameters:
        * `Whether to return the orientations in degrees or radians
        """
        
        # Store as radians/degrees and return as radians/degrees
        if (self.degrees and degrees) or (not self.degrees and not degrees):
            return self.phi_1, self.Phi, self.phi_2
        
        # Store as radians but return as degrees
        elif not self.degrees and degrees:
            return rad_to_deg(self.phi_1), rad_to_deg(self.Phi), rad_to_deg(self.phi_2)
        
        # Store as degrees but return as radians
        elif self.degrees and not degrees:
            return deg_to_rad(self.phi_1), deg_to_rad(self.Phi), deg_to_rad(self.phi_2)

    def set_grain_id(self, grain_id:int) -> int:
        """
        Sets the grain ID
        
        Parameters:
        * `grain_id`: The ID of the grain to which the element belongs
        """
        self.grain_id = grain_id
        
    def get_grain_id(self) -> int:
        """
        Returns the grain ID
        """
        return self.grain_id

def get_void_element(grain_id:int=VOID_ELEMENT_ID) -> Element:
    """
    Returns a void element
    
    Parameters:
    * `grain_id`: The ID of the element
    """
    return Element(*VOID_ORIENTATION, grain_id=grain_id)

def initialise_plot(mesh, figure_x:float, exodus_face:str="xy", custom_size:tuple=None, custom_limits:tuple=None) -> None:
    """
    Initialises the plot

    Parameters:
    * `mesh`:        The mesh object
    * `figure_x`:    Factor to change the size of the figure
    * `exodus_face`: The face to be plotted
    """
    
    # Get mesh bounds
    x_min, x_max = get_exodus_bounds(mesh, exodus_face[0])
    y_min, y_max = get_exodus_bounds(mesh, exodus_face[1])

    # Initialises the figure
    if custom_size != None:
        plt.figure(figsize=custom_size)
    else:
        plt.figure(figsize=(figure_x, (y_max-y_min)/(x_max-x_min)*figure_x))

    # Defines the limits
    if custom_limits != None:
        plt.xlim(custom_limits[0])
        plt.ylim(custom_limits[1])
    else:
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
    
    # Apply additional formatting
    plt.gca().set_aspect("equal", adjustable="datalim")
    plt.gca().invert_yaxis()

def plot_mesh(mesh, grain_dict:dict, element_indexes:list, ipf:str="x",
              directions:str="xy", positive:bool=True) -> None:
    """
    Creates a plot of the mesh

    Parameters:
    * `mesh`:            The mesh object
    * `exodus_path`:     Path to the exodus file
    * `grain_dict`:      Map of exodus elements to the initial EBSD map
    * `element_indexes`: Indexes of elements to plot
    * `thickness`:       Number of elements along the z-axis
    * `ipf`:             The IPF direction to plot the mesh
    * `directions`:      The directions of to plot the mesh
    * `positive`:        Whether to plot the positive or negative face
    """

    # Initialise perspective
    dir_to_index = {"x": 0, "y": 1, "z": 2}
    included = [dir_to_index[dir] for dir in directions]
    excluded = (set([0,1,2])-set(included)).pop()
    _, max_value = get_exodus_bounds(mesh, "xyz"[excluded])
    face = 0 if positive else max_value

    # Iterate and plot each grain
    for i, grain in enumerate(mesh):
        for cell_id in element_indexes[i]:

            # Get grain ID
            exo_id = int(str(mesh.get_block_name(i)).split(" ")[-1])
            if not exo_id in grain_dict:
                continue
            
            # Get cell coordinates and ignore if not surface
            cell_coordinates = get_surface_cell_coordinates(grain, cell_id, included, excluded, face)
            if cell_coordinates == []:
                continue

            # Get IPF colour
            element = grain_dict[exo_id][cell_id]
            orientation = element.get_orientation(degrees=True)
            colour = [rgb/255 for rgb in euler_to_rgb(*orientation, ipf=ipf)]
            
            # Plot the cell
            cell_coordinates = order_vertices(cell_coordinates)
            if cell_coordinates == None:
                continue
            # polygon = patches.Polygon(cell_coordinates, closed=True, fill=True, facecolor=colour, edgecolor="black")
            polygon = patches.Polygon(cell_coordinates, closed=True, fill=True, facecolor=colour, edgecolor=None)
            plt.gca().add_patch(polygon)

def get_exodus_bounds(mesh, direction:str="z") -> float:
    """
    Gets the dimension of the mesh in a certain direction

    Parameters:
    * `mesh`:      The mesh object
    * `direction`: The direction ("x", "y", "z")

    Returns the dimension
    """
    bounds = mesh.bounds # (xmin, xmax, ymin, ymax, zmin, zmax)
    direction_map = {"x": (0, 1), "y": (2, 3), "z": (4, 5)}
    min_bound, max_bound = direction_map[direction]
    return bounds[min_bound], bounds[max_bound]

def get_element_indexes(mesh, directions:str, positive:bool) -> list:
    """
    Gets the indexes of the elements in a mesh that lie on a face

    Parameters:
    * `mesh`:        The mesh object
    * `directions`:  The directions of to plot the mesh
    * `positive`:    Whether to plot the positive or negative face
    * `exodus_path`: Path to the exodus file

    Returns the list of coordinates for the surface of the mesh
    """

    # Get information about the mesh
    dir_to_index = {"x": 0, "y": 1, "z": 2}
    included = [dir_to_index[dir] for dir in directions]
    excluded = (set([0,1,2])-set(included)).pop()
    _, max_value = get_exodus_bounds(mesh, "xyz"[excluded])
    face = 0 if positive else max_value

    # Iterate through grains
    element_index_grid = []
    for i in range(mesh.n_blocks):
        grain = mesh[i]
        element_index_list = []

        # Iterate through elements in grains
        for cell_id in range(grain.n_cells):
            point_ids = grain.get_cell(cell_id).point_ids
            coordinates = [list(point) for point in grain.points[point_ids]]
            for coordinate in coordinates:
                if coordinate[excluded] == face:
                    element_index_list.append(cell_id)
                    break
        element_index_grid.append(element_index_list)
    
    # Return element indexes
    return element_index_grid

def get_surface_cell_coordinates(grain:pv.core.pointset.UnstructuredGrid, cell_id:int,
                                 included:list, excluded:int, face:bool) -> list:
    """
    Gets the coordinates of a cell

    Parameters:
    * `grain`:    The grain object
    * `cell_id`:  The cell ID
    * `included`: The coordinate positions to plot
    * `excluded`: The coordinate position to not plot
    * `face`:     Whether to plot the positive or negative face

    Returns the list of coordinates for the surface corners of the cell
    """
    point_ids = grain.get_cell(cell_id).point_ids
    coordinates = [list(point) for point in grain.points[point_ids]]
    coordinates = sorted(coordinates, key=lambda c: c[excluded], reverse=face)[:4] # doesn't work for major deformation
    coordinates = [[c[i] for i in included] for c in coordinates]
    return coordinates

def on_segment(p:list, q:list, r:list): 
    """
    Check if point q lies on the line segment pr;
    only for colinear points

    Parameters:
    * `p`: The first point
    * `q`: The second point
    * `r`: The three point
    """
    return (
        (q[0] <= max(p[0], r[0])) and
        (q[0] >= min(p[0], r[0])) and
        (q[1] <= max(p[1], r[1])) and
        (q[1] >= min(p[1], r[1]))
    )
  
def orientation(p:list, q:list, r:list): 
    """
    Gets the orientations of three points

    Parameters:
    * `p`: The first point
    * `q`: The second point
    * `r`: The three point
    """
    val = (float(q[1] - p[1]) * (r[0] - q[0])) - (float(q[0] - p[0]) * (r[1] - q[1])) 
    if (val > 0): 
        return 1
    elif (val < 0): 
        return 2
    else: 
        return 0
  
def do_intersect(p_1:list, q_1:list, p_2:list, q_2:list) -> bool: 
    """
    Checks whether two segments intersect

    Parameters:
    * `p_1`: The first point of the first segment
    * `q_1`: The second point of the first segment
    * `p_2`: The first point of the second segment
    * `q_2`: The second point of the second segment

    Returns whether an intersection occurs or not
    """
    o_1 = orientation(p_1, q_1, p_2) 
    o_2 = orientation(p_1, q_1, q_2) 
    o_3 = orientation(p_2, q_2, p_1) 
    o_4 = orientation(p_2, q_2, q_1) 
    return (
        ((o_1 != o_2) and (o_3 != o_4)) or
        ((o_1 == 0) and on_segment(p_1, p_2, q_1)) or
        ((o_2 == 0) and on_segment(p_1, q_2, q_1)) or
        ((o_3 == 0) and on_segment(p_2, p_1, q_2)) or
        ((o_4 == 0) and on_segment(p_2, q_1, q_2))
    )

def order_vertices(unordered_vertices:list) -> list:
    """
    Orders a list of four vertices so that
    they form a quadrilateral

    Parameters:
    * `unordered_vertices`: List of unordered vertices

    Returns the ordered lists of vertices
    """
    for permutation in itertools.permutations(unordered_vertices[1:]):
        vertices = [unordered_vertices[0]] + list(permutation)
        intersect = do_intersect(vertices[0], vertices[2], vertices[1], vertices[3])
        if intersect:
            return vertices

# Calls the main function
if __name__ == "__main__":
    main()

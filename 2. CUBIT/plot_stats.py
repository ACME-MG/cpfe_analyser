"""
 Title:         Statistics Plotter
 Description:   Plots the statistics of meshes
 Author:        Janzen Choi

"""

# Libraries
from exodus import get_exodus_grain_areas

# Constants
MESH_PATH = "data/617_s3/10u/mesh.e"

print(get_exodus_grain_areas(MESH_PATH))

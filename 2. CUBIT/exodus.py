"""
 Title:         Exodus
 Description:   For manipulating exodus files
 Author:        Janzen Choi

"""

# Libraries
import math
import netCDF4 as nc
import pyvista as pv

def get_exodus_length(exodus_path:str, direction:str="z") -> float:
    """
    Gets the dimension of the mesh in a certain direction

    Parameters:
    * `exodus_path`: Path to the exodus file
    * `direction`:   The direction to scale ("x", "y", "z")

    Returns the dimension
    """
    ds = nc.Dataset(exodus_path, mode="r")
    coordinates = ds.variables[f"coord{direction}"][:]
    ds.close()
    dimension = max(coordinates) - min(coordinates)
    return dimension

def get_equiv_radii(exo_path:str) -> list:
    """
    Calculates the equivalent radius of the grains of a mesh
    
    Parameters:
    * `exo_path`: Path to the exodus file
    
    Returns the radii as a list
    """
    
    # Calculate the areas
    thickness = get_exodus_length(exo_path, "z")
    area_list = []
    mesh = pv.read(exo_path)[0]
    for grain in mesh:
        area_list.append(grain.volume / thickness)
    
    # Calculate the equivalent radii and return
    radius_list = [math.pow(area/math.pi, 0.5) for area in area_list]
    return radius_list

def get_circularity(exo_path:str) -> list:
    """
    Calculates the circularity of grains in a mesh
    
    Parameters:
    * `exo_path`: Path to the exodus file
    
    Returns the circularities as a list
    """
    
    # Initialise
    mesh = pv.read(exo_path)[0]
    thickness = get_exodus_length(exo_path, "z")
    circularity_list = []
    
    # Iterate through grains
    for grain in mesh:
        
        # Get all edges in the grain
        edge_list = []
        for cell_id in range(grain.n_cells):
            cell = grain.get_cell(cell_id)
            for edge in cell.edges:
                points = [list(point) for point in edge.points]
                if points[0][-1] != 0 or points[1][-1] != 0:
                    continue
                points = tuple(tuple(point[:-1]) for point in points)
                edge_list.append(points)
            
        # Identify boundary cells
        edge_count = {}
        for edge in edge_list:
            sorted_edge = tuple(sorted(edge))
            if sorted_edge in edge_count:
                edge_count[sorted_edge] += 1
            else:
                edge_count[sorted_edge] = 1
        boundary_edges = [edge for edge, count in edge_count.items() if count == 1]
        
        # Calculate perimeter
        distances   = [math.sqrt(math.pow(be[0][0]-be[1][0], 2) + math.pow(be[0][1]-be[1][1], 2)) for be in boundary_edges]
        perimeter   = sum(distances)
        area        = grain.volume / thickness
        circularity = 4*math.pi*area / perimeter**2
        circularity_list.append(circularity)
        
    # Return
    return circularity_list
        
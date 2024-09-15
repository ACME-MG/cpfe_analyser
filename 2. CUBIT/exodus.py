"""
 Title:         Exodus
 Description:   For manipulating exodus files
 Author:        Janzen Choi

"""

# Libraries
import netCDF4 as nc

def get_exodus_grain_areas(exodus_path:str) -> list:
    """
    Calculates the grain areas of a mesh
    
    Parameters:
    * `exodus_path`: Path to the exodus file
    
    Returns the grain areas as a list
    """
    return [1,2,3]

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

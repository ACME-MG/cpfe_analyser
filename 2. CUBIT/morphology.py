"""
 Title:         Morphology
 Description:   Analyses the morphology of different meshes
 Author:        Janzen Choi

"""

# Libraries
from exodus import get_equiv_radii, get_circularity
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import lognorm

def round_sf(value:float, sf:int) -> float:
    """
    Rounds a float to a number of significant figures

    Parameters:
    * `value`: The value to be rounded; accounts for lists
    * `sf`:    The number of significant figures

    Returns the rounded number
    """
    if isinstance(value, list):
        return [round_sf(v, sf) for v in value]
    format_str = "{:." + str(sf) + "g}"
    rounded_value = float(format_str.format(value))
    return rounded_value

def initialise_plot(x_label:str="", y_label:str="", x_max:float=None, y_max:float=None) -> None:
    """
    Initialises a plot
    
    Parameters:
    * `x_max`:   Upper limit for the x-axis
    * `y_max`:   Upper limit for the y-axis
    * `x_label`: Label for the x-axis
    * `y_label`: Label for the y-axis
    """
    plt.figure(figsize=(5,5))
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":")
    plt.xlabel(x_label, fontsize=11) if x_label != "" else None
    plt.ylabel(y_label, fontsize=11) if y_label != "" else None
    plt.xlim(0,x_max) if x_max != None else None
    plt.ylim(0,y_max) if y_max != None else None
    
def format_and_save_plot(plot_path:str) -> None:
    """
    Formats and saves a plot
    
    Parameters:
    * `plot_path`: Path to save the plot
    """
    plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white")
    plt.savefig(plot_path)
    plt.cla()
    plt.clf()
    plt.close()

def plot_distribution(value_list:list, colour:str, label:str, min_value:float=None,
                      max_value:float=None, settings:dict={}) -> None:
    """
    Plots a distribution plot given a list of values
    
    Parameters:
    * `value_list`: List of values
    * `colour`:     Colour to use when plotting
    * `label`:      Label to use when plotting
    * `min_value`:  Minimum value to plot the distribution
    * `max_value`:  Maximum value to plot the distribution
    """
    shape, loc, scale = lognorm.fit(value_list)
    min_value = min(value_list) if min_value == None else min_value
    max_value = max(value_list) if max_value == None else max_value
    x_list = np.linspace(min_value, max_value, 1000)
    y_list = lognorm.pdf(x_list, shape, loc, scale)
    mu = round_sf(np.log(scale), 3)
    mu = str(mu) + "0" if len(str(mu)) <= 3 else mu
    sigma = round_sf(shape, 3)
    plt.plot(x_list, y_list, colour, label=f"{label} ({mu}, {sigma})", **settings)

# Define mesh information
mesh_info_list = [
    {"path": "data/617_s3_z1/5um/mesh.e",  "colour": "silver",  "label": "5µm  "},
    {"path": "data/617_s3_z1/10um/mesh.e", "colour": "orange",  "label": "10µm"},
    {"path": "data/617_s3_z1/15um/mesh.e", "colour": "gold",    "label": "15µm"},
    {"path": "data/617_s3_z1/20um/mesh.e", "colour": "brown",   "label": "20µm"},
    {"path": "data/617_s3_z1/25um/mesh.e", "colour": "red",     "label": "25µm"},
    {"path": "data/617_s3_z1/30um/mesh.e", "colour": "magenta", "label": "30µm"},
    {"path": "data/617_s3_z1/35um/mesh.e", "colour": "purple",  "label": "35µm"},
    {"path": "data/617_s3_z1/40um/mesh.e", "colour": "blue",    "label": "40µm"},
    {"path": "data/617_s3_z1/45um/mesh.e", "colour": "cyan",    "label": "45µm"},
    {"path": "data/617_s3_z1/50um/mesh.e", "colour": "green",   "label": "50µm"},
]

# Plot equivalent radius distribution
initialise_plot("Equivalent Radius", "Probability Density", 250, 0.05)
for i, mesh_info in enumerate(mesh_info_list):
    settings = {"linewidth": 3} if i == 0 else {"linestyle": "dashed"}
    radius_list = get_equiv_radii(mesh_info["path"])
    radius_list = sorted(radius_list, reverse=True)[2:] # exclude grips
    plot_distribution(radius_list, mesh_info["colour"], mesh_info["label"], 0, 250, settings)
format_and_save_plot("eq_rad.png")

# Plot circularity distribution
initialise_plot("Circularity", "Probability Density", 3.0, 7.0)
for i, mesh_info in enumerate(mesh_info_list):
    settings = {"linewidth": 3} if i == 0 else {"linestyle": "dashed"}
    circularity_list = get_circularity(mesh_info["path"])
    circularity_list = sorted(circularity_list, reverse=True)[2:] # exclude grips
    plot_distribution(circularity_list, mesh_info["colour"], mesh_info["label"], 0.0, 3.0, settings)
format_and_save_plot("circularity.png")

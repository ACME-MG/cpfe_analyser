"""
 Title:         Number of Grains
 Description:   Plots the number of grains across multiple meshes
 Author:        Janzen Choi

"""

# Libraries
import matplotlib.pyplot as plt
import sys; sys.path += [".."]
from __common__.plotter import save_plot
import matplotlib.colors as mcolors

# Mesh information
MESH_INFO_LIST = [
    {"path": "data/617_s3_z1/5um",  "colour": "silver",     "label": "5µm  "},
    {"path": "data/617_s3_z1/10um", "colour": "tab:orange", "label": "10µm"},
    {"path": "data/617_s3_z1/15um", "colour": "gold",       "label": "15µm"},
    {"path": "data/617_s3_z1/20um", "colour": "sienna",     "label": "20µm"},
    {"path": "data/617_s3_z1/25um", "colour": "tab:red",    "label": "25µm"},
    {"path": "data/617_s3_z1/30um", "colour": "magenta",    "label": "30µm"},
    {"path": "data/617_s3_z1/35um", "colour": "tab:purple", "label": "35µm"},
    {"path": "data/617_s3_z1/40um", "colour": "tab:blue",   "label": "40µm"},
    {"path": "data/617_s3_z1/45um", "colour": "tab:green",  "label": "45µm"},
]

# Other constants
TICK_SIZE  = 12
LABEL_SIZE = 14
INCREMENT  = 5
WIDTH      = 4

def main():
    """
    Main function
    """

    # Calculate resolutions
    resolution_list = [mesh_info["label"] for mesh_info in MESH_INFO_LIST]
    resolution_list = [int(r.replace("µm","")) for r in resolution_list]

    # Calculate the number of grains in each mesh
    num_grains_list = []
    for mesh_info in MESH_INFO_LIST:
        with open(f"{mesh_info['path']}/grain_map.csv", "r") as fh:
            num_grains = len(fh.readlines())-2
            num_grains_list.append(num_grains)

    # Plot number of grains
    plt.figure(figsize=(5,5), dpi=200)
    plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
    plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
    plt.xlabel("Resolution (µm)", fontsize=LABEL_SIZE)
    plt.ylabel("Number of grains", fontsize=LABEL_SIZE)
    for resolution, num_grains, colour in zip(resolution_list, num_grains_list, [mi["colour"] for mi in MESH_INFO_LIST]):
        colour = lighten_color(colour, 0.5)
        plt.bar([resolution], [num_grains], color=colour, zorder=3, width=WIDTH, edgecolor="black")
    padding = INCREMENT-WIDTH/2
    plt.xlim(min(resolution_list)-padding, max(resolution_list)+padding)
    plt.ylim(0,600)
    print(num_grains_list)
    plt.xticks(resolution_list, fontsize=TICK_SIZE)
    plt.yticks(fontsize=TICK_SIZE)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    save_plot("results/num_grains.png")

def lighten_color(colour:str, amount=0.5):
    """
    Lightens a MatplotLib colour

    Parameters:
    * `colour`: The name of the colour
    * `amount`: The amount to lighten
    """
    c = mcolors.to_rgba(colour)  # Convert to RGBA
    white = (1, 1, 1, 1)  # White color
    return tuple((1 - amount) * c[i] + amount * white[i] for i in range(4))

# Main function
if __name__ == "__main__":
    main()

"""
 Title:         Number of Grains
 Description:   Plots the number of grains across multiple meshes
 Author:        Janzen Choi

"""

# Libraries
import matplotlib.pyplot as plt
import numpy as np
import sys; sys.path += [".."]
from __common__.io import csv_to_dict
from __common__.plotter import Plotter, save_plot

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
COLOUR     = "black"

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
    plt.xlim(min(resolution_list)-2.5, max(resolution_list)+2.5)
    plt.ylim(0,600)
    plt.plot(resolution_list, num_grains_list, color=COLOUR, linewidth=3, marker="s")
    print(num_grains_list)
    plt.xticks(resolution_list, fontsize=TICK_SIZE)
    plt.yticks(fontsize=TICK_SIZE)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    save_plot("results/num_grains.png")

# Main function
if __name__ == "__main__":
    main()

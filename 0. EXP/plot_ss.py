"""
 Title:        Plot Experimental
 Description:  Generates plots of the experimental data
 Author:       Janzen Choi

"""

# Libraries
# import matplotlib.pyplot as plt
import sys; sys.path += [".."]
import matplotlib.pyplot as plt
import math, numpy as np
from __common__.io import  read_excel
from __common__.general import remove_nan
from __common__.plotter import Plotter, save_plot

# Paths
ENG_SS_PATH  = f"/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/data/2024-06-26 (ansto_617_s3)/j3_sscurve.xlsx"
ENG_SS_SHEET = "curve"
ENG_SI_SHEET = "ebsd_points"
TRU_SS_PATH  = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/data/2024-06-26 (ansto_617_s3)/sscurve_corrected_janzen_3.xlsx"
TRU_SS_SHEET = "Sheet1"
RAW_RT_PATH  = "data/617_s3_20um_raw_rt.csv"
PRC_RT_PATH  = "data/617_s3_20um_prc_rt.csv"

# Colours
ENG_COLOUR    = "silver"
TRU_COLOUR    = "gray"
UNLOAD_COLOUR = "tab:blue"

def main():
    """
    Main function
    """

    # Read engineering stress-strain data
    strain_list = [0] + read_excel(ENG_SS_PATH, ENG_SS_SHEET, 0)[1:]
    stress_list = [0] + read_excel(ENG_SS_PATH, ENG_SS_SHEET, 1)[1:]
    keep_list = [i for i in range(len(strain_list)) if isinstance(strain_list[i], float)]
    strain_list = [strain_list[i] for i in keep_list]
    stress_list = [stress_list[i] for i in keep_list]

    # Read unloading stress-strain data
    unload_strain_list = remove_nan(read_excel(ENG_SS_PATH, ENG_SI_SHEET, 1))
    unload_stress_list = remove_nan(read_excel(ENG_SS_PATH, ENG_SI_SHEET, 2))
    
    # Calculate unloaded stress-strain points
    youngs_index = 10
    youngs = stress_list[youngs_index]/strain_list[youngs_index]
    unloaded_strain_list, unloaded_stress_list = get_unloaded(unload_strain_list, unload_stress_list, youngs, 15)
    
    # Add unloading sections into stress-strain data
    for i in range(len(unload_strain_list)):
        strain_list += list(np.linspace(unload_strain_list[i], unloaded_strain_list[i], 50))
        stress_list += list(np.linspace(unload_stress_list[i], unloaded_stress_list[i], 50))

    # Read true stress-strain data
    true_strain_list = [0] + read_excel(TRU_SS_PATH, TRU_SS_SHEET, 5)[1:]
    true_stress_list = [0] + read_excel(TRU_SS_PATH, TRU_SS_SHEET, 6)[1:]
    true_unload_strain_list = remove_nan(read_excel(TRU_SS_PATH, TRU_SS_SHEET, 14)[1:])
    true_unload_stress_list = remove_nan(read_excel(TRU_SS_PATH, TRU_SS_SHEET, 15)[1:])

    # Plot stress-strain curve
    plotter = Plotter("Strain", "Stress", "mm/mm", "MPa")
    plotter.prep_plot(size=14)
    plotter.set_limits((0,0.7), (0,1600))
    plt.scatter(strain_list,             stress_list,             color=ENG_COLOUR, s=8**2)
    plt.plot   (true_strain_list,        true_stress_list,        color=TRU_COLOUR, linewidth=3)
    plt.scatter(unloaded_strain_list,    unloaded_stress_list,    color=UNLOAD_COLOUR, s=4**2, zorder=3, marker="s")
    plt.scatter(true_unload_strain_list, true_unload_stress_list, color="None", s=8**2, linewidths=2, zorder=3, marker="o", edgecolor=UNLOAD_COLOUR)
    handles = [
        plt.scatter([], [], color=ENG_COLOUR,    label="Eng. SS", s=8**2),
        plt.plot   ([], [], color=TRU_COLOUR,    label="True SS", linewidth=3)[0],
        plt.scatter([], [], color=UNLOAD_COLOUR, label="Eng. EBSD", marker="s", s=4**2),
        plt.scatter([], [], color="None",        label="True EBSD", s=8**2, linewidths=2, marker="o", edgecolor=UNLOAD_COLOUR)
    ]
    legend = plt.legend(handles=handles, framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
    plt.gca().add_artist(legend)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(1)
    save_plot("results/exp_ss.png")

def truify(strain_list:list, stress_list:list) -> tuple:
    """
    Converts engineering strain and stress into true strain and stress

    Parameters:
    * `strain_list`: Engineering strains
    * `stress_list`: Engineering stresses

    Returns true strain and stress as lists
    """
    true_strain_list, true_stress_list = [], []
    for i in range(len(strain_list)):
        true_strain = math.log(1+strain_list[i])
        true_stress = stress_list[i]*(1+strain_list[i])
        true_strain_list.append(true_strain)
        true_stress_list.append(true_stress)
    return true_strain_list, true_stress_list

def get_unloaded(unload_strain_list:list, unload_stress_list:list, youngs:float, unloaded_stress:float) -> tuple:
    """
    Simulates unloading sections based on the unloading intervals

    Parameters:
    * `unload_strain_list`: List of strain values where the specimen was unloaded
    * `unload_stress_list`: List of stress values where the specimen was unloaded
    * `youngs`:             Youngs modulus
    * `unloaded_stress`:    Stress to unload to
    
    Returns the unloaded strain and stress values as a tuple
    """
    unloaded_strain_list = []
    for unload_strain, unload_stress in zip(unload_strain_list, unload_stress_list):
        unloaded_strain = unload_strain-(unload_stress-unloaded_stress)/youngs
        unloaded_strain_list.append(unloaded_strain)
    return unloaded_strain_list, [unloaded_stress]*len(unloaded_strain_list)

# Main functionn caller
if __name__ == "__main__":
    main()

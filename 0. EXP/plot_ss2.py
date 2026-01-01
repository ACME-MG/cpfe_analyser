"""
 Title:        Plot Experimental
 Description:  Generates plots of the experimental data
 Author:       Janzen Choi

"""

# Libraries
# import matplotlib.pyplot as plt
import sys; sys.path += [".."]
import matplotlib.pyplot as plt
import math
from __common__.io import  read_excel, csv_to_dict
from __common__.general import remove_nan, get_thinned_list
from __common__.plotter import Plotter, save_plot

# Paths
TRU_SS_PATH  = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/data/2024-06-26 (ansto_617_s3)/sscurve_corrected_janzen_3.xlsx"
TRU_SS_SHEET = "Sheet1"

# Colours
INL_COLOUR    = "silver"
TRU_COLOUR    = "gray"
UNLOAD_COLOUR = "tab:red"

def main():
    """
    Main function
    """

    # Get INL tensile data
    inl_dict = csv_to_dict("data/AirBase_20_D5.csv")
    inl_strain_list = get_thinned_list(inl_dict["strain"], 500)[:-1]
    inl_stress_list = get_thinned_list(inl_dict["stress"], 500)[:-1]

    # Read true stress-strain data
    tru_dict = csv_to_dict("data/617_s3_40um_exp.csv")
    tru_strain_list = tru_dict["strain"]
    tru_stress_list = tru_dict["stress"]
    unl_strain_list = remove_nan(read_excel(TRU_SS_PATH, TRU_SS_SHEET, 14)[1:])
    unl_stress_list = remove_nan(read_excel(TRU_SS_PATH, TRU_SS_SHEET, 15)[1:])

    # Plot stress-strain curve
    plotter = Plotter("Strain", "Stress", "mm/mm", "MPa")
    plotter.prep_plot(size=14)
    plotter.set_limits((0,0.5), (0,1800))
    plt.scatter(inl_strain_list, inl_stress_list, color=INL_COLOUR, s=8**2)
    plt.scatter(tru_strain_list, tru_stress_list, color=TRU_COLOUR, s=6**2)
    plt.scatter(unl_strain_list, unl_stress_list, color="None", linewidths=2, marker="o", edgecolor=UNLOAD_COLOUR, s=8**2, zorder=3)
    handles = [
        plt.scatter([], [], color=INL_COLOUR,    label="Standard Specimen", s=8**2),
        plt.scatter([], [], color=TRU_COLOUR,    label="Miniaturised Specimen", s=6**2),
        plt.scatter([], [], color="None", label="EBSD Data Acquisition", linewidths=2, marker="o", edgecolor=UNLOAD_COLOUR, s=8**2),
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
    true_strain_list, tru_stress_list = [], []
    for i in range(len(strain_list)):
        true_strain = math.log(1+strain_list[i])
        true_stress = stress_list[i]*(1+strain_list[i])
        true_strain_list.append(true_strain)
        tru_stress_list.append(true_stress)
    return true_strain_list, tru_stress_list

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

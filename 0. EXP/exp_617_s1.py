"""
 Title:        Experimental Grouper
 Description:  For grouping reorientation data with stress strain data
 Author:       Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import numpy as np
import matplotlib.pyplot as plt
from __common__.general import round_sf
from __common__.io import csv_to_dict, dict_to_csv
from __common__.interpolator import Interpolator

# Paths
TENSILE_PATH       = "./data/617_s1_tc.csv"
REORIENTATION_PATH = "./data/617_s1_reorientation.csv"
EXPERIMENTAL_PATH  = "./results/617_s1_exp.csv"

# Fields
TIME_FIELD   = "time"   # s
STRESS_FIELD = "stress" # MPa
STRAIN_FIELD = "strain" # mm/mm

# Other constants
NUM_POINTS = 100
RAW_STRAINS = [0.0, 0.003, 0.008, 0.022, 0.033, 0.036, 0.048, 0.075, 0.094, 0.110, 0.132, 0.144, 0.165, 0.189]
NEXT_STRAIN = 0.20

def process_tensile(strain_list:list, stress_list:list, num_points:int=40) -> list:
    """
    For converting in-situ EBSD tensile data into normal tensile data

    Parameters:
    * `strain_list`: The list of strain values
    * `stress_list`: The list of stress values
    * `num_points`:  Number of points to reduce the lists to
    
    Returns the lists of indexes corresponding to the chosen strain and stress values
    """

    # Initialise
    interval_size = len(strain_list) // num_points
    index_list = []

    # Get indexes
    for i in range(1,num_points+1):
        lower_index = max(0, interval_size*i-interval_size//2)
        upper_index = min(len(strain_list)-1, interval_size*i+interval_size//2)
        interval_indexes = range(lower_index, upper_index+1)
        stress_interval = [stress_list[j] for j in interval_indexes]
        stress_index = stress_interval.index(max(stress_interval))
        index_list.append(stress_index + min(interval_indexes))

    # Return
    return index_list

# Read files
tensile_dict = csv_to_dict(TENSILE_PATH)
reorient_dict = csv_to_dict(REORIENTATION_PATH)

# Convert tensile data
index_list   = process_tensile(tensile_dict[STRAIN_FIELD], tensile_dict[STRESS_FIELD])
strain_list  = [tensile_dict[STRAIN_FIELD][i] for i in index_list]
stress_list  = [tensile_dict[STRESS_FIELD][i] for i in index_list]
interpolator = Interpolator(strain_list, stress_list)
strain_list  = [round_sf(strain, 5) for strain in RAW_STRAINS + list(np.linspace(NEXT_STRAIN, max(strain_list), NUM_POINTS-len(RAW_STRAINS)))]
stress_list  = [0] + [round_sf(stress, 5) for stress in interpolator.evaluate(strain_list[1:])]
time_list    = [round_sf(time, 5) for time in np.linspace(0, max(tensile_dict[TIME_FIELD]), NUM_POINTS)]

# Extend reorientation data
for key in reorient_dict.keys():
    reorient_dict[key] = reorient_dict[key] + [reorient_dict[key][-1]]*(len(strain_list)-len(reorient_dict[key]))

# Create tensile dictionary
tc_dict = {
    "time":   time_list,
    "strain": strain_list,
    "stress": stress_list,
}

# Create dictionary of other information
other_dict = {
    "temperature": 20,     # C
    "strain_rate": 1e-4,   # /s
    "youngs":      211000, # MPa
    "poissons":    0.30,
    "type":        "tensile",
    "title":       "617_s1",
}

# Combine dictionaries and save
exp_dict = {**tc_dict, **other_dict, **reorient_dict}
dict_to_csv(exp_dict, EXPERIMENTAL_PATH)

# Plot tensile curve
plt.figure(figsize=(5, 5))
plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
plt.scatter(tensile_dict[STRAIN_FIELD], tensile_dict[STRESS_FIELD], color="grey", label="Raw"),
plt.plot(strain_list, stress_list, color="red",  label="Processed"),
legend = plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
plt.gca().add_artist(legend)
plt.xlim(0,1.20)
plt.ylim(0,1200)
plt.savefig(f"results/617_s1_ss.png")


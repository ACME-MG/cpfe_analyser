"""
 Title:        Experimental Grouper
 Description:  For grouping reorientation data with stress strain data
 Author:       Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import matplotlib.pyplot as plt
from __common__.general import round_sf, get_thinned_list, remove_nan, get_closest
from __common__.io import csv_to_dict, dict_to_csv, read_excel

# Input Paths
SS_PATH  = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/data/2024-06-26 (ansto_617_s3)/sscurve_corrected_janzen_3.xlsx"
# PHI_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/ebsd_mapper/240813093129_617_s3_mapping/reorientation.csv"
# PHI_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/ebsd_mapper/240912165226_617_s3_10u/reorientation.csv"
PHI_PATH = "data/617_s3_20um_reorientation.csv"

# Output Paths
EXPERIMENTAL_PATH = "./results/617_s3_exp.csv"

# Other Constants
THIN_AMOUNT = 500

# Read stress-strain information
time_list   = [0] + read_excel(SS_PATH, "Sheet1", 0)[1:] # exclude header
strain_list = [0] + read_excel(SS_PATH, "Sheet1", 5)[1:] # exclude header
stress_list = [0] + read_excel(SS_PATH, "Sheet1", 6)[1:] # exclude header

# Read stress-strain at EBSD
strain_intervals = read_excel(SS_PATH, "Sheet1", 14)[1:] # exclude header
strain_intervals = remove_nan(strain_intervals)
stress_intervals = read_excel(SS_PATH, "Sheet1", 15)[1:] # exclude header
stress_intervals = remove_nan(stress_intervals)
time_intervals   = [get_closest(strain_list, time_list, strain) for strain in strain_intervals]

# Thin the stress-strain information
time_list   = get_thinned_list(time_list, THIN_AMOUNT)
strain_list = get_thinned_list(strain_list, THIN_AMOUNT)
stress_list = get_thinned_list(stress_list, THIN_AMOUNT)

# Read and extend reorientation data
phi_dict = csv_to_dict(PHI_PATH)
for phi in phi_dict.keys():
    phi_dict[phi] = [round_sf(np, 5) for np in phi_dict[phi]]
    phi_dict[phi] = phi_dict[phi][1:]

# Create new stress-strain dictionary
new_ss_dict = {
    "time":   [round_sf(time, 5)   for time in time_list],
    "strain": [round_sf(strain, 5) for strain in strain_list],
    "stress": [round_sf(stress, 5) for stress in stress_list],
}

# Create dictionary for EBSD intervals
interval_dict = {
    "time_intervals":   [round_sf(ti, 5) for ti in time_intervals],
    "strain_intervals": [round_sf(si, 5) for si in strain_intervals],
    "stress_intervals": [round_sf(si, 5) for si in stress_intervals],
}

# Create dictionary of other information
other_dict = {
    "temperature": 20,     # C
    "strain_rate": 1e-4,   # /s
    "youngs":      211000, # MPa
    "poissons":    0.30,
    "type":        "tensile",
    "title":       "617_s3",
}

# Combine dictionaries and save
exp_dict = {**new_ss_dict, **interval_dict, **phi_dict, **other_dict}
dict_to_csv(exp_dict, EXPERIMENTAL_PATH)

# Plot tensile curve
plt.figure(figsize=(5, 5))
plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
plt.scatter(strain_list, stress_list, color="grey", label="Raw"),
plt.xlim(0,0.5)
plt.ylim(0,1400)
plt.savefig(f"results/617_s3_ss.png")

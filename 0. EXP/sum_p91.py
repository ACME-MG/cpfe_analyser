"""
 Title:        Summariser
 Description:  Gets necessary information about deformation behaviour and summarises
 Author:       Janzen Choi

"""

# Libraries
import math
import sys; sys.path += [".."]
from __common__.general import round_sf
from __common__.io import csv_to_dict, dict_to_csv, read_excel
from __common__.orientation import deg_to_rad

# Index to control which data to access
SAMPLE_INDEX = int(sys.argv[1])

# Main directory containing data
DATA_DIR = "/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/data/20240516 (ondrej_P91)"

# Input for tensile
TENSILE_PATH   = f"{DATA_DIR}/Bulk vs Micro tensile plots.xlsx"
TENSILE_SHEET  = "Micro Tensile_P91"
TENSILE_STRAIN = [0,4,8][SAMPLE_INDEX]
TENSILE_STRESS = [1,5,9][SAMPLE_INDEX]
TENSILE_TIME   = [2,6,10][SAMPLE_INDEX]

# Input for grain files
FILE_NAME = "grainsExportColumnsTable.csv"
GRAIN_START_PATH = [
    f"{DATA_DIR}/S1_2/P91_Eurofer97P91UNIRS1_2MapData22/{FILE_NAME}",
    f"{DATA_DIR}/S2_2/P91_Eurofer97P91UNIRS2_2MapData23/{FILE_NAME}",
    f"{DATA_DIR}/S3_2/P91_Eurofer97P91UNIRS3_2MapData24/{FILE_NAME}",
][SAMPLE_INDEX]
GRAIN_END_PATH = [
    f"{DATA_DIR}/S1_2/P91_Eurofer97P91UNIRS1_postMapData53/{FILE_NAME}",
    f"{DATA_DIR}/S2_2/P91_Eurofer97P91UNIRS1_postMapData53/{FILE_NAME}",
    f"{DATA_DIR}/S3_2/P91_Eurofer97P91UNIRS3_postMapData55/{FILE_NAME}",
][SAMPLE_INDEX]

# Input for grain information
GRAIN_PHI_1  = "Euler_phi1"
GRAIN_PHI    = "Euler_Phi"
GRAIN_PHI_2  = "Euler_phi2"
GRAIN_WEIGHT = "grainNumPixels"

# Define grain map
GRAIN_MAP_0 = {1: 7, 2: 1, 3: 9, 4: 10, 5: 5, 6: 6, 7: 8, 10: 21, 13: 11, 14: 15, 16: 35, 17: 12, 18: 33,
               19: 25, 20: 23, 21: 37, 22: 24, 23: 29, 24: 38, 25: 39, 26: 40, 28: 42, 30: 47, 31: 48, 32: 51,
               33: 55, 36: 41, 37: 44, 39: 45, 40: 52, 41: 53, 42: 61, 44: 46, 45: 49, 48: 54, 50: 60, 51: 66}
GRAIN_MAP_1 = {}
GRAIN_MAP_2 = {1: 4, 2: 5, 4: 7, 3: 6, 5: 8, 6: 14, 7: 9, 9: 15, 10: 16, 11: 13, 12: 17, 13: 20, 14: 11, 15: 10,
               16: 22, 19: 30, 20: 27, 21: 29, 22: 32, 23: 19, 24: 31, 29: 26, 30: 37, 31: 28, 32: 35, 33: 33,
               35: 38, 36: 41, 37: 40, 38: 47, 40: 34, 41: 44, 43: 52, 44: 51, 45: 58, 47: 39, 48: 42, 49: 46,
               50: 43, 52: 53, 53: 54, 54: 55, 55: 57, 57: 59, 65: 60, 68: 69, 70: 71, 71: 72, 72: 75, 77: 77,
               79: 80, 80: 81, 81: 82, 82: 83}
GRAIN_MAP = [GRAIN_MAP_0, GRAIN_MAP_1, GRAIN_MAP_2][SAMPLE_INDEX]

# Output Parameters
OUTPUT_FOLDER  = "results"
OUTPUT_MAPPING = f"{OUTPUT_FOLDER}/p91_s{SAMPLE_INDEX+1}_map"
OUTPUT_TENSILE = f"{OUTPUT_FOLDER}/p91_s{SAMPLE_INDEX+1}_exp"
OUTPUT_GRAINS  = f"{OUTPUT_FOLDER}/p91_s{SAMPLE_INDEX+1}_grains"

# Other Constants
YOUNGS   = 190000 # MPa
POISSONS = 0.28

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

def get_trajectory(strain_list:list, start:float, end:float) -> list:
    """
    Gets the linear trajectory an orientation

    Parameters:
    * `strain_list`: Strains
    * `start`:      Orientation starting value (deg)
    * `end`:        Orientation ending value (deg)
    
    Returns the orientation values (rad)
    """
    min_strain = min(strain_list)
    max_strain = max(strain_list)
    gradient = (end-start)/(max_strain-min_strain)
    linear = lambda x : gradient*(x-min_strain) + start
    orientation_list = [linear(strain) for strain in strain_list]
    orientation_list = deg_to_rad(orientation_list)
    return orientation_list

# Read and process tensile data
time_list   = read_excel(TENSILE_PATH, TENSILE_SHEET, TENSILE_TIME)
strain_list = read_excel(TENSILE_PATH, TENSILE_SHEET, TENSILE_STRAIN)
stress_list = read_excel(TENSILE_PATH, TENSILE_SHEET, TENSILE_STRESS)
strain_list, stress_list = truify(strain_list, stress_list)
strain_rate = round_sf(max(strain_list)/max(time_list), 7)

# Read and process start CSV
grain_start_dict   = csv_to_dict(GRAIN_START_PATH)
grain_start_phi_1  = grain_start_dict[GRAIN_PHI_1]
grain_start_Phi    = grain_start_dict[GRAIN_PHI]
grain_start_phi_2  = grain_start_dict[GRAIN_PHI_2]
grain_start_weight = grain_start_dict[GRAIN_WEIGHT]

# Read and process end CSV
grain_end_dict   = csv_to_dict(GRAIN_END_PATH)
grain_end_phi_1  = grain_end_dict[GRAIN_PHI_1]
grain_end_Phi    = grain_end_dict[GRAIN_PHI]
grain_end_phi_2  = grain_end_dict[GRAIN_PHI_2]
grain_end_weight = grain_end_dict[GRAIN_WEIGHT]

# Get grain trajectories of grains we can map
tensile_grain_dict = {}
for start_index in GRAIN_MAP.keys():
    end_index = GRAIN_MAP[start_index]
    tensile_grain_dict[f"g{start_index}_phi_1"] = get_trajectory(strain_list, grain_start_phi_1[start_index-1], grain_end_phi_1[end_index-1])
    tensile_grain_dict[f"g{start_index}_Phi"]   = get_trajectory(strain_list, grain_start_Phi[start_index-1],   grain_end_Phi[end_index-1])
    tensile_grain_dict[f"g{start_index}_phi_2"] = get_trajectory(strain_list, grain_start_phi_2[start_index-1], grain_end_phi_2[end_index-1])

# Create dictionaries for tensile data
tensile_curve_dict = {
    "time":  time_list,
    "strain": strain_list,
    "stress": stress_list,
}
tensile_field_dict = {
    "temperature": 20,
    "type":       "tensile",
    "medium":     "Air",
    "strain_rate": strain_rate,
    "youngs":     YOUNGS,
    "poissons":   POISSONS,
}

# Package and save tensile data
tensile_dict = {**tensile_curve_dict, **tensile_grain_dict, **tensile_field_dict}
dict_to_csv(tensile_dict, f"{OUTPUT_TENSILE}.csv")

# Package and save grain data
grain_dict = {
    "phi_1": deg_to_rad(grain_start_phi_1),
    "Phi":   deg_to_rad(grain_start_Phi),
    "phi_2": deg_to_rad(grain_start_phi_2),
    "weight": grain_start_weight
}
dict_to_csv(grain_dict, f"{OUTPUT_GRAINS}.csv", False)

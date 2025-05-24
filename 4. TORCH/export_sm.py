"""
 Title:         Plot Surrogate
 Description:   Plots the response of the surrogate model
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += ["..", "/home/janzen/code/mms"]
from __common__.io import dict_to_csv
from __common__.surrogate import Model

# Paths
EXP_PATH = "data/617_s3_exp.csv"
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo"
MOOSE_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/moose_sim"
DIRECTORY = f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)/250318012854_i21_surrogate"
SUR_PATH = f"{DIRECTORY}/sm.pt"
MAP_PATH = f"{DIRECTORY}/map.csv"
RESULTS_PATH = "../5. ASMBO/data/summary.csv"

# Constants
MAX_STRAIN = 0.29
PARAM_LIST = [1938.4, 0.17367, 149.03, 6.1855]

# Main function
def main():
    model = Model(SUR_PATH, MAP_PATH, EXP_PATH, MAX_STRAIN)
    res_dict = model.get_response(PARAM_LIST)
    for key in res_dict.keys():
        res_dict[key] = list(res_dict[key])
    res_dict["average_strain"] = res_dict["strain"]
    res_dict["average_stress"] = res_dict["stress"]
    dict_to_csv(res_dict, RESULTS_PATH)

# Calls the main function
if __name__ == "__main__":
    main()

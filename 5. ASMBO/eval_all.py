"""
 Title:         Evaluates a bunch of simulations together
 Description:   Plots the errors for a set of simulations with Pareto efficiency and dominance
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
from pareto import pareto
from tolerance import tolerance

# Evaluation function
EVALUATE = pareto
# EVALUATE = tolerance

# Paths
ASMBO_DIR_DICT = {
    2: [
        "2025-01-23 (vh_sm2_i32)",
    ],
    4: [
        "2025-01-20 (vh_sm4_i17)",
        "2025-01-20 (vh_sm4_i22)",
        "2025-01-20 (vh_sm4_i25)",
        "2025-01-25 (vh_sm4_i29)",
    ],
    6: [
        "2025-01-22 (vh_sm6_i23)",
        "2025-01-22 (vh_sm6_i18)",
    ],
    8: [
        "2025-01-18 (vh_sm8_i24)",
        "2025-01-19 (vh_sm8_i22)",
        "2025-01-25 (vh_sm8_i16)",
    ],
    10: [
        "2025-01-24 (vh_sm10_i22)",
        "2025-01-26 (vh_sm10_i10)",
    ],
    12: [
        "2025-01-26 (vh_sm12_i15)",
    ],
    14: [
        "2025-01-26 (vh_sm14_i25)",
        "2025-01-27 (vh_sm14_i23)",
    ],
    16: [
        "2025-01-21 (vh_sm16_i19)",
    ],
}
SIM_DATA_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo"
EXP_DATA_PATH = "data/617_s3_40um_exp.csv"
RESULTS_PATH  = "results"

# Fields
STRAIN_FIELD = "average_strain"
STRESS_FIELD = "average_stress"

# Grain information
CAL_GRAIN_IDS = [59, 63, 86, 237, 303]
VAL_GRAIN_IDS = [44, 53, 60, 78, 190]

# Error information
ERROR_COLOUR = "black"
PD_COLOUR    = "tab:red"
PE_COLOUR    = "tab:green"
KP_COLOUR    = "tab:blue"

def paretos():
    """
    Main function
    """
    for init_evals in ASMBO_DIR_DICT.keys():
        sim_path_list = [f"{SIM_DATA_PATH}/{sim_dir}" for sim_dir in ASMBO_DIR_DICT[init_evals]]
        num_evals = [EVALUATE(sim_path) for sim_path in sim_path_list]
        print(init_evals, num_evals)

# Calls the main function
if __name__ == "__main__":
    paretos()

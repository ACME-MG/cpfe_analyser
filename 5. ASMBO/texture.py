"""
 Title:         Plot Texture
 Description:   Plots the texture of the optimised simulations at certain strains
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
from __common__.general import transpose
from __common__.io import csv_to_dict, dict_to_csv
from __common__.interpolator import intervaluate

# Paths
EXP_PATH     = "data/617_s3_40um_exp.csv"
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo"
MOOSE_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/moose_sim"
# RESULTS_PATH = "results/pf"
RESULTS_PATH = "/mnt/c/Users/janzen/Desktop/texture/scripts/data/temp/"

# Simulation information
SIM_INFO_LIST = [
    
    # VH 
    # f"{ASMBO_PATH}/2025-03-09 (vh_pin2_sm8_i25)/250308143546_i4_simulate",
    # f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310145710_i22_simulate",
    # f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)/250318014435_i21_simulate",
    # f"{ASMBO_PATH}/2025-03-25 (vh_x_sm8_i31)/250325072901_i16_simulate",
    # f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310161708_i25_simulate",
    # ===================================================================
    # f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)/250318014435_i21_simulate",
    # f"{MOOSE_PATH}/2025-03-15 (617_s3_vh_x_hr)",

    # # LH2
    # f"{ASMBO_PATH}/2025-03-25 (lh2_x_sm8_i19)/250323214745_i7_simulate",
    # f"{ASMBO_PATH}/2025-03-28 (lh2_x_sm8_i29)/250327093649_i16_simulate",
    # f"{ASMBO_PATH}/2025-03-31 (lh2_x_sm8_i31)/250330063457_i18_simulate",
    # f"{ASMBO_PATH}/2025-04-02 (lh2_x_sm8_i23)/250401051359_i10_simulate",
    # f"{ASMBO_PATH}/2025-03-31 (lh2_x_sm8_i31)/250330213453_i29_simulate",
    # ===================================================================
    # f"{ASMBO_PATH}/2025-03-28 (lh2_x_sm8_i29)/250327093649_i16_simulate",
    # f"{MOOSE_PATH}/2025-04-05 (617_s3_lh2_di_x_hr)",

    # # LH6
    f"{ASMBO_PATH}/2025-04-09 (lh6_x_sm8_i44)/250407052902_i6_simulate",
    f"{ASMBO_PATH}/2025-04-14 (lh6_x_sm8_i32)/250413031321_i23_simulate",
    f"{ASMBO_PATH}/2025-04-18 (lh6_x_sm8_i27)/250418123844_i27_simulate",
    f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250422034348_i36_simulate",
    f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250420224600_i20_simulate",
    # ===================================================================
    # f"{ASMBO_PATH}/2025-04-09 (lh6_x_sm8_i44)/250407052902_i6_simulate",
    # f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250422034348_i36_simulate",
    # f"{MOOSE_PATH}/2025-04-28 (617_s3_lh6_di_x_hr)",

    # VH CV1
    # f"{ASMBO_PATH}/2025-06-03 (vh_x_sm8_i97_val)/250601032309_i20_simulate",
    # f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i13_val)/250609093847_i7_simulate",
    # # f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310145710_i22_simulate",
    # # f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i52_val)/250608004520_i6_simulate",
    # f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i29_cv2)/250707105726_i9_simulate",
    # f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i28_val)/250608021510_i6_simulate",
    # f"{ASMBO_PATH}/2025-06-11 (vh_x_sm8_i52_val)/250609222901_i8_simulate",
    
    # VH CV2
    # f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i14_cv2)/250708043252_i4_simulate",
    # f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i14_cv2)/250708020621_i1_simulate",
    # f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i26_cv2)/250707055737_i4_simulate",
    # f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i29_cv2)/250707093811_i6_simulate",
    # f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i30_cv2)/250708071533_i13_simulate",
    # f"{ASMBO_PATH}/2025-07-08 (vh_x_sm8_i32_cv2)/250708043852_i9_simulate",
]

# Plotting parameters
# EVAL_STRAINS = [0.01, 0.052, 0.106, 0.208, 0.29]
EVAL_STRAINS = [0.0001, 0.00063414, 0.00153, 0.00494, 0.0098, 0.01483, 0.02085, 0.02646, 0.03516, 0.04409, 0.05197, 0.06013, 0.07059, 0.08208, 0.09406, 0.10561, 0.11929, 0.13656, 0.15442, 0.18237, 0.20849, 0.23627, 0.26264, 0.28965]
FIELD_LIST   = ["phi_1", "Phi", "phi_2", "volume", "area"]

def main():
    """
    Main function
    """

    # Initialise
    get_grain_ids    = lambda dict : [int(key.replace("g","").replace("_phi_1","")) for key in dict.keys() if "_phi_1" in key]
    get_trajectories = lambda dict, grain_ids : [transpose([dict[f"g{grain_id}_{phi}"] for phi in FIELD_LIST if f"g{grain_id}_{phi}" in dict.keys()]) for grain_id in grain_ids]

    # Read experimental data
    exp_dict = csv_to_dict(EXP_PATH)
    exp_strain_list = exp_dict["strain_intervals"]
    exp_grain_ids = get_grain_ids(exp_dict)
    exp_trajectories = get_trajectories(exp_dict, exp_grain_ids)

    # Read simulation data
    sim_dict_list = [csv_to_dict(f"{si}/summary.csv") for si in SIM_INFO_LIST]
    for sim_dict in sim_dict_list:
        sim_grain_ids = get_grain_ids(sim_dict)
        sim_dict["trajectories"] = get_trajectories(sim_dict, sim_grain_ids)
    
    # Extract orientations
    exp_orientations_list = []
    sim_orientations_list_list = []
    for i, eval_strain in enumerate(EVAL_STRAINS):
        exp_orientations = [intervaluate_orientation(exp_strain_list, et, eval_strain) for et in exp_trajectories]
        sim_orientations_list = []
        for sim_dict in sim_dict_list:
            sim_orientations = [intervaluate_orientation(sim_dict["average_strain"], st, eval_strain) for st in sim_dict["trajectories"]]
            sim_orientations_list.append(sim_orientations)
        exp_orientations_list.append(exp_orientations)
        sim_orientations_list_list.append(sim_orientations_list)

    # Export data
    for i in range(len(EVAL_STRAINS)):
        exp_orientations = exp_orientations_list[i]
        sim_orientations_list = sim_orientations_list_list[i]
        export_data(exp_orientations, f"exp_i{i+1}")
        for j, sim_orientations in enumerate(sim_orientations_list):
            export_data(sim_orientations, f"sim_s{j+1}_i{i+1}")

def export_data(euler_list:list, file_name:str) -> None:
    """
    Exports texture data as CSV files

    Parameters:
    * `euler_list`: The list of orientations in euler-bunge form (rads)
    * `file_name`:  Name of file to save the data (without extension)
    """
    phi_lists = transpose(euler_list)
    phi_dict = {}
    for i, field in enumerate(FIELD_LIST):
        if i < len(phi_lists):
            phi_dict[field] = phi_lists[i]
    dict_to_csv(phi_dict, f"{RESULTS_PATH}/{file_name}.csv")

def intervaluate_orientation(strain_list:list, orientation_list:list, strain:float) -> list:
    """
    Interpolates and evaluates from a list of orientations

    Parameters:
    * `strain_list`:      List of strain values
    * `orientation_list`: List of orientation components
    * `strain`:           Strain value to evaluate

    Returns the intervaluated orientation
    """
    components = transpose(orientation_list)
    new_components = [intervaluate(strain_list, component, strain) for component in components]
    return transpose(new_components)

# Main function
if __name__ == "__main__":
    main()

"""
 Title:         Plot Optimised Simulations
 Description:   Plots the response of the optimised simulations
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import numpy as np
from __common__.io import csv_to_dict
from __common__.general import round_sf, transpose
from __common__.interpolator import intervaluate
from __common__.orientation import euler_to_quat, get_geodesic, rad_to_deg
from scipy.spatial.transform import Rotation

# Constants
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/H0419460/results/asmbo"
MODEL_INDEX = int(sys.argv[1]) if len(sys.argv) > 1 else 0

# VH information
if MODEL_INDEX == 0:
    SIM_PATHS = [
        f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310145710_i22_simulate",
        f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310145710_i22_simulate",
        f"{ASMBO_PATH}/2025-03-18 (vh_x_sm8_i41)/250318014435_i21_simulate",
        f"{ASMBO_PATH}/2025-03-25 (vh_x_sm8_i31)/250325072901_i16_simulate",
        f"{ASMBO_PATH}/2025-03-10 (vh_pin2_sm8_i25)/250310161708_i25_simulate",
    ]
    SIM_TI = "data/vh_ti.csv"

# LH2 information
if MODEL_INDEX == 1:
    SIM_PATHS = [
        f"{ASMBO_PATH}/2025-03-25 (lh2_x_sm8_i19)/250323214745_i7_simulate",
        f"{ASMBO_PATH}/2025-03-28 (lh2_x_sm8_i29)/250327093649_i16_simulate",
        f"{ASMBO_PATH}/2025-03-31 (lh2_x_sm8_i31)/250330063457_i18_simulate",
        f"{ASMBO_PATH}/2025-04-02 (lh2_x_sm8_i23)/250401051359_i10_simulate",
        f"{ASMBO_PATH}/2025-03-31 (lh2_x_sm8_i31)/250330213453_i29_simulate",
    ]
    SIM_TI = "data/lh2_ti.csv"

# LH6 information
if MODEL_INDEX == 2:
    SIM_PATHS = [
        f"{ASMBO_PATH}/2025-04-09 (lh6_x_sm8_i44)/250407052902_i6_simulate",
        f"{ASMBO_PATH}/2025-04-14 (lh6_x_sm8_i32)/250413031321_i23_simulate",
        f"{ASMBO_PATH}/2025-04-18 (lh6_x_sm8_i27)/250418123844_i27_simulate",
        f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250422034348_i36_simulate",
        f"{ASMBO_PATH}/2025-04-23 (lh6_x_sm8_i51)/250420224600_i20_simulate",
    ]
    SIM_TI = "data/lh6_ti.csv"

# Grain IDs
GRAIN_IDS = [
    [14, 72, 95, 101, 207, 240, 262, 287],  # Calibration
    [39, 50, 138, 164, 185, 223, 243, 238], # Validation
]

# Main function
def main():

    # Prepare data
    exp_dict = csv_to_dict("data/617_s3_40um_exp.csv")
    eval_strains = exp_dict["strain_intervals"][1:]
    sim_dict_list = [csv_to_dict(f"{sp}/summary.csv") for sp in SIM_PATHS]
    sim_ti_dict = csv_to_dict(SIM_TI)

    # euler_list = get_euler_list(sim_dict_list, GRAIN_IDS[0][0], eval_strains[-1])
    # rt_sd = get_euler_sd(euler_list)
    # print(rad_to_deg(rt_sd))

    # Get SD of stress-strain data
    stress_grid = [[intervaluate(sim_dict["average_strain"], sim_dict["average_stress"], eval_strain)
                    for sim_dict in sim_dict_list] for eval_strain in eval_strains]
    stress_sd = get_avg_sd(stress_grid)
    print(f"SS:  {round_sf(stress_sd,3)} MPa")

    # Get SD of reorientation trajectories
    for i, grain_ids in enumerate(GRAIN_IDS):
        rt_sd = get_euler_avg_sd(sim_dict_list, grain_ids, eval_strains)
        print(f"RT{i+1}: {round_sf(rad_to_deg(rt_sd),3)}°")

    # Get SD of texture evolution
    all_grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in sim_dict_list[0].keys() if "_phi_1" in key]
    te_sd = get_euler_avg_sd(sim_dict_list, all_grain_ids, eval_strains)
    print(f"TE:  {round_sf(rad_to_deg(te_sd),3)}°")

    # Get SD of texture indexes
    ti_grid = transpose([sim_ti_dict[f"Run {i+1}"] for i in range(5)])[1:]
    ti_sd = get_avg_sd(ti_grid)
    print(f"TI:  {round_sf(ti_sd,3)}")
    
    # Separator
    print("="*50)

def get_euler_avg_sd(sim_dict_list:list, grain_ids:list, eval_strains:list) -> float:
    """
    Gets the average standard deviation of orientations

    Parameters:
    * `sim_dict_list`: List of simulation information
    * `grain_ids`:     List of grain IDs to evaluate
    * `eval_strains`:  List of evaluated strains

    Returns the euler angles at a strain
    """
    euler_sd_list = []
    for eval_strain in eval_strains:
        for grain_id in grain_ids:
            euler_list = get_euler_list(sim_dict_list, grain_id, eval_strain)
            euler_sd = get_euler_sd(euler_list)
            euler_sd_list.append(euler_sd)
    euler_avg_sd = np.average(euler_sd_list)
    return euler_avg_sd

def get_euler_sd(euler_list:list, axes:str="zxz") -> float:
    """
    Calculates the standard deviation given a list of euler angles (rads)

    Parameters:
    * `euler_list`: The list of euler angles

    Returns the standard deviation
    """
    
    # Get average euler angle (Fréchet)
    rotations = Rotation.from_euler(axes, euler_list, degrees=False)
    euler_avg = np.mod(rotations.mean().as_euler(axes, degrees=False), 2*np.pi)
    
    # Calculate geodesic distances
    quat_avg = euler_to_quat(euler_avg)
    geodesic_list = []
    for euler in euler_list:
        quat = euler_to_quat(euler)
        geodesic = get_geodesic(quat, quat_avg)
        geodesic_list.append(geodesic**2)
    
    # Calculate standard deviation
    geodesic_sd = np.sqrt(np.sum(geodesic_list)/(len(geodesic_list)-1))
    return geodesic_sd

def get_euler_list(sim_dict_list:list, grain_id:int, eval_strain:float) -> list:
    """
    Gets the orientation information at a given strain

    Parameters:
    * `sim_dict_list`: List of simulation information
    * `grain_id`:      The grain ID to evaluate
    * `eval_strain`:   The evaluated strain

    Returns the euler angles at a strain
    """
    euler_list = []
    for sim_dict in sim_dict_list:
        strain_list = sim_dict["average_strain"]
        phi_list = [sim_dict[f"g{grain_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]
        euler = [intervaluate(strain_list, phi, eval_strain) for phi in phi_list]
        euler_list.append(euler)
    return euler_list

def get_avg_sd(data_grid:list) -> float:
    """
    Gets the average standard deviation

    Parameters:
    * `data_grid`: A list of lists of data

    Returns the average standard deviation
    """
    data_sd_list = []
    for data_list in data_grid:
        data_var = np.var(data_list, ddof=1)
        data_sd = np.sqrt(data_var)
        data_sd_list.append(data_sd)
    return np.average(data_sd_list)

# Calls the main function
if __name__ == "__main__":
    main()

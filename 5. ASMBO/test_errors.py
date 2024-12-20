"""
 Title:         Test Errors
 Description:   Compares the errors used in the optimisation and analysis
 Author:        Janzen Choi

"""

# Libraries
import math, numpy as np
import sys; sys.path += [".."]
from __common__.io import csv_to_dict
from __common__.interpolator import Interpolator
from __common__.orientation import get_geodesic, euler_to_quat
from __common__.analyse import get_geodesics, get_stress

# Constants
GRAIN_IDS   = [59, 63, 86, 237, 303]
EXP_PATH    = "data/617_s3_40um_exp.csv"
SIM_DIR     = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim"
SIM_PATH    = f"{SIM_DIR}/2024-11-05 (617_s3_40um_lh2_opt)/summary.csv"
NUM_POINTS  = 50
MAX_STRAIN  = 0.1
EVAL_X_LIST = list(np.linspace(0, MAX_STRAIN, NUM_POINTS))

def main():
    """
    Main function
    """

    # Read experimental and simulated data
    exp_data = csv_to_dict(EXP_PATH)
    sim_data = csv_to_dict(SIM_PATH)

    # Compare stress errors
    print("Stress error")
    opt_stress = get_opt_stress(exp_data, sim_data)
    anl_stress = get_stress(exp_data["stress"], sim_data["average_stress"], exp_data["strain"], sim_data["average_strain"], EVAL_X_LIST)
    print(opt_stress, anl_stress)

    # Compare geodesic errors
    print("Geodesic errors")
    for grain_id in GRAIN_IDS:
        opt_geodesic = get_opt_geodesic(exp_data, sim_data, grain_id)
        geodesic_grid = get_geodesics([grain_id], exp_data, sim_data, exp_data["strain_intervals"], sim_data["average_strain"], EVAL_X_LIST)
        anl_geodesic = np.average([np.sqrt(np.average([g**2 for g in gg])) for gg in geodesic_grid])
        print(opt_geodesic, anl_geodesic)

def get_opt_stress(exp_data:dict, prd_data:dict) -> float:
    """
    Gets the stress error using the optimiser's implementationn

    Parameters:
    * `exp_data`: Experimental data
    * `prd_data`: Predicted data

    Returns the stress error
    """
    
    # Get experimental data
    x_list           = exp_data["strain"]
    y_list           = exp_data["stress"]
    exp_interpolator = Interpolator(x_list, y_list, len(x_list))
    exp_x_end        = min(x_list[-1], MAX_STRAIN) if MAX_STRAIN != None else x_list[-1]
    exp_y_end        = exp_interpolator.evaluate([exp_x_end])[0]
    avg_abs_y        = np.average([abs(y) for y in y_list if y < exp_y_end])

    # Get predicted data
    x_label    = "average_strain"
    y_label    = "average_stress"
    try:
        prd_interpolator = Interpolator(prd_data[x_label], prd_data[y_label], len(prd_data[x_label]))
    except:
        return
    
    # Calculate error
    prd_x_list = np.linspace(prd_data[x_label][0], min(exp_x_end, prd_data[x_label][-1]), NUM_POINTS)
    prd_y_list = prd_interpolator.evaluate(prd_x_list)
    exp_y_list = exp_interpolator.evaluate(prd_x_list)
    try:
        area = [math.pow(prd_y_list[i] - exp_y_list[i], 2) for i in range(NUM_POINTS) if prd_x_list[i] <= exp_x_end]
        return math.sqrt(np.average(area)) / avg_abs_y
    except OverflowError:
        return

def get_opt_geodesic(exp_data:dict, prd_data:dict, grain_id:int) -> float:
    """
    Gets the geodesic error using the optimiser's implementationn

    Parameters:
    * `exp_data`: Experimental data
    * `prd_data`: Predicted data
    * `grain_id`: ID of grain

    Returns the geodesic error
    """

    # Get strain list
    labels = ["strain_intervals", f"g{grain_id}_phi_1", f"g{grain_id}_Phi", f"g{grain_id}_phi_2"]
    exp_x_list = exp_data[labels[0]]

    # Interpolate experimental orientations
    exp_itp_list = []
    for label in labels[1:]:
        exp_phi_list = exp_data[label]
        exp_itp = Interpolator(exp_x_list, exp_phi_list, len(exp_phi_list))
        exp_itp_list.append(exp_itp)
    
    # Calculate experimental orientations at target strain values
    exp_phi_list = []
    for eval_x in EVAL_X_LIST:
        exp_phi = [exp_itp.evaluate([eval_x])[0] for exp_itp in exp_itp_list]
        exp_phi_list.append(exp_phi)

    # Get predicted data
    prd_x_list = prd_data["average_strain"]
    
    # Interpolate predicted orientations
    prd_itp_list = []
    for label in labels[1:]:
        prd_phi_list = prd_data[label]
        prd_itp = Interpolator(prd_x_list, prd_phi_list, len(prd_phi_list))
        prd_itp_list.append(prd_itp)
        
    # Calculate predicted orientations at target strain values
    prd_phi_list = []
    for eval_x in EVAL_X_LIST:
        prd_phi = [prd_itp.evaluate([eval_x])[0] for prd_itp in prd_itp_list]
        prd_phi_list.append(prd_phi)
    
    # Calculate geodesic distances
    geodesic_list = []
    for exp_phi, prd_phi in zip(exp_phi_list, prd_phi_list):
        exp_quat = euler_to_quat(exp_phi)
        prd_quat = euler_to_quat(prd_phi)
        geodesic = get_geodesic(exp_quat, prd_quat)
        geodesic_list.append(geodesic)
    
    # Average and return
    square_geodesic = [math.pow(geodesic, 2) for geodesic in geodesic_list]
    return math.sqrt(np.average(square_geodesic))

# Main function caller
if __name__ == "__main__":
    main()

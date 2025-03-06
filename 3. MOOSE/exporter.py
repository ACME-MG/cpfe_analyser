"""
 Title:         Exporter
 Description:   Creates a CSV file of the texture at a specific strain
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
from __common__.general import round_sf
from __common__.interpolator import intervaluate
from __common__.io import csv_to_dict, dict_to_csv

# Paths
ASMBO_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo"
MOOSE_PATH = "/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/moose_sim"
DATA_PATH = f"{MOOSE_PATH}/2025/2025-02-09 (617_s3_10um_lh2)"

# Evaluation information
STRAIN_FIELD = "average_strain"
EVAL_STRAIN = 0.30

def main():
    """
    Main function
    """

    # Read data
    data_dict = csv_to_dict(f"{DATA_PATH}/summary.csv")
    strain_list = data_dict[STRAIN_FIELD]
    grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in data_dict.keys() if "_phi_1" in key]

    # Intervaluate the orientations
    phi_1_list = [intervaluate(strain_list, data_dict[f"g{grain_id}_phi_1"], EVAL_STRAIN) for grain_id in grain_ids]
    Phi_list   = [intervaluate(strain_list, data_dict[f"g{grain_id}_Phi"], EVAL_STRAIN) for grain_id in grain_ids]
    phi_2_list = [intervaluate(strain_list, data_dict[f"g{grain_id}_phi_2"], EVAL_STRAIN) for grain_id in grain_ids]

    # Get the volumes
    volume_list = [intervaluate(strain_list, data_dict[f"g{grain_id}_volume"], EVAL_STRAIN) for grain_id in grain_ids]

    # Create new dictionary
    new_dict = {
        "grain_id": grain_ids,
        "phi_1":    round_sf(phi_1_list, 5),
        "Phi":      round_sf(Phi_list, 5),
        "phi_2":    round_sf(phi_2_list, 5),
        "volume":   round_sf(volume_list, 5),
    }

    # Save as CSV
    dict_to_csv(new_dict, "results/texture.csv")

# Calls the main function
if __name__ == "__main__":
    main()

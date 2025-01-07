"""
 Title:         ID mapper
 Description:   Maps the grain IDs based on previous grain ID maps
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
from __common__.io import csv_to_dict

# Constants
GRAIN_MAP     = "data/res_grain_map.csv"
SRC_EBSD_ID   = "ebsd_4" # 20um
DST_EBSD_ID   = "ebsd_2" # 10um
SRC_GRAIN_IDS = [59, 63, 86, 237, 303]

def main():
    """
    Main function
    """

    # Read grain IDs
    grain_map = csv_to_dict(GRAIN_MAP)
    ebsd_src_ids = grain_map[SRC_EBSD_ID]
    ebsd_dst_ids = grain_map[DST_EBSD_ID]
    
    # Use grain mapping and print
    dst_grain_ids = [ebsd_dst_ids[ebsd_src_ids.index(grain_id)] for grain_id  in SRC_GRAIN_IDS]
    dst_grain_ids = [int(grain_id) for grain_id in dst_grain_ids]
    print(dst_grain_ids)

# Main function caller
if __name__ == "__main__":
    main()

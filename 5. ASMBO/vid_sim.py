"""
 Title:         Video Simulations
 Description:   Summarises the responses for a set of simulations through video
 Author:        Janzen Choi

"""

# Libraries
import cv2, os
import matplotlib.pyplot as plt
import numpy as np
import sys; sys.path += [".."]
from __common__.general import transpose
from __common__.io import csv_to_dict
from __common__.pole_figure import get_lattice, IPF
from __common__.plotter import define_legend

# Constants
ASMBO_DIR     = "2025-01-09 (lh_0p3_i16)"
SIM_DATA_PATH = f"/mnt/c/Users/janzen/OneDrive - UNSW/PhD/results/asmbo/{ASMBO_DIR}"
EXP_DATA_PATH = "data/617_s3_40um_exp.csv"
RESULTS_PATH  = "results"
STRAIN_FIELD  = "average_strain"
STRESS_FIELD  = "average_stress"
# GRAIN_IDS     = [59, 63, 86, 237, 303]
GRAIN_IDS     = [44, 56, 60, 141, 207]
EXP_COLOUR    = "silver"
SIM_COLOUR    = "tab:red"
SIM_LABEL     = "Calibration"

def main():
    """
    Main function
    """

    # Gets experimental data
    exp_dict = csv_to_dict(EXP_DATA_PATH)

    # Get summary of simulations
    dir_path_list = [f"{SIM_DATA_PATH}/{dir_path}" for dir_path in os.listdir(SIM_DATA_PATH) if "sim" in dir_path]
    sum_path_list = [f"{dir_path}/summary.csv" for dir_path in dir_path_list if os.path.exists(f"{dir_path}/summary.csv")]
    sim_dict_list = [csv_to_dict(summary_path) for summary_path in sum_path_list]

    # Plot and record stress-strain response
    vw_ss = get_video_writer(f"{RESULTS_PATH}/sum_ss.mp4", (1000,1000))
    for sim_dict in sim_dict_list:
        plt.figure(figsize=(5,5), dpi=200)
        plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
        plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":", alpha=0.5)
        plt.xlabel("Strain (mm/mm)")
        plt.ylabel("Stress (MPa)")
        plt.xlim(0,0.5)
        plt.ylim(0,1400)
        plt.scatter(exp_dict["strain"], exp_dict["stress"], s=8**2, color=EXP_COLOUR, label="Experiment")
        plt.plot(sim_dict[STRAIN_FIELD], sim_dict[STRESS_FIELD], linewidth=3, color=SIM_COLOUR, label=SIM_LABEL)
        plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="lower right")
        write_to_video(vw_ss)
    vw_ss.release()

    # Plot and record reorientation trajectories
    vw_rt = get_video_writer(f"{RESULTS_PATH}/sum_rt.mp4", (1200,1000))
    direction = [1,0,0]
    get_trajectories = lambda dict, grain_ids : [transpose([dict[f"g{grain_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]) for grain_id in grain_ids]
    for sim_dict in sim_dict_list:
        ipf = IPF(get_lattice("fcc"))
        exp_trajectories = get_trajectories(exp_dict, GRAIN_IDS)
        ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": EXP_COLOUR, "linewidth": 3})
        ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": EXP_COLOUR, "head_width": 0.01, "head_length": 0.015})
        ipf.plot_ipf_trajectory([[et[0]] for et in exp_trajectories], direction, "scatter", {"color": EXP_COLOUR, "s": 8**2})
        for exp_trajectory, grain_id in zip(exp_trajectories, GRAIN_IDS):
            ipf.plot_ipf_trajectory([[exp_trajectory[0]]], direction, "text", {"color": "black", "fontsize": 8, "s": grain_id})
        cal_trajectories = get_trajectories(sim_dict, GRAIN_IDS)
        ipf.plot_ipf_trajectory(cal_trajectories, direction, "plot", {"color": SIM_COLOUR, "linewidth": 2, "zorder": 3})
        ipf.plot_ipf_trajectory(cal_trajectories, direction, "arrow", {"color": SIM_COLOUR, "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
        ipf.plot_ipf_trajectory([[ct[0]] for ct in cal_trajectories], direction, "scatter", {"color": SIM_COLOUR, "s": 6**2, "zorder": 3})
        define_legend([EXP_COLOUR, SIM_COLOUR], ["Experiment", SIM_LABEL], ["scatter", "line"], fontsize=12)
        write_to_video(vw_rt)
    vw_rt.release()

def get_video_writer(path:str, size:tuple) -> cv2.VideoWriter:
    """
    Gets a video writer
    
    Parameters:
    * `path`: Path to save the video
    * `size`: Size of the video frames

    Returns a video writer object
    """
    frame_rate = 1
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    video_writer = cv2.VideoWriter(path, fourcc, frame_rate, size)
    return video_writer

def write_to_video(video_writer:cv2.VideoWriter):
    """
    Gets the current figure and adds it to a video writer
    
    Parameters:
    * `video_writer`: The video writer
    """
    figure = plt.gcf()
    figure.canvas.draw()
    image = np.frombuffer(figure.canvas.tostring_rgb(), dtype=np.uint8)
    image = image.reshape(figure.canvas.get_width_height()[::-1] + (3,))
    image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)
    video_writer.write(image)
    plt.close(figure)

# Calls the main function
if __name__ == "__main__":
    main()


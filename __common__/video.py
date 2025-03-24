"""
 Title:         Video
 Description:   Converts plots into videos
 Author:        Janzen Choi

"""

# Libraries
import cv2
import numpy as np
from matplotlib.figure import Figure

# Constants

# The Video class
class Video:

    def __init__(self, video_name:str, frame_rate:float, frame_size:tuple):
        """
        Constructor for the Video class;
        creates a video writer object
        
        Parameters:
        * `video_name`: Name to use to save the video file (without extension)
        * `frame_rate`: Number of frames per second
        * `frame_size`: Dimensions of the video
        """
        self.frame_rate = frame_rate
        self.frame_size = frame_size
        fourcc = cv2.VideoWriter_fourcc(*"mp4v") # codec for mp4
        self.video_writer = cv2.VideoWriter(f"{video_name}.mp4", fourcc, frame_rate, frame_size)

    def write_to_video(self, figure:Figure) -> None:
        """
        Adds a frame to the video

        Parameters:
        * `figure`: The figure object
        """
        figure.canvas.draw()
        image = np.frombuffer(figure.canvas.tostring_rgb(), dtype=np.uint8)
        image = image.reshape(figure.canvas.get_width_height()[::-1] + (3,))
        image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)
        self.video_writer.write(image)

    def __del__(self) -> None:
        """
        Destructor for the Video class;
        releases the video writer object
        """
        self.video_writer.release()


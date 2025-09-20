from scipy.spatial.transform import Rotation as R
import numpy as np

# sample Euler angles (ZYX convention)
angles = np.random.randn(100, 3) * 0.1  

rots = R.from_euler('ZYX', angles)

# Convert to rotation vectors (log map)
rotvecs = rots.as_rotvec()  # shape (N,3)

# covariance in tangent space
cov = np.cov(rotvecs.T)
print(cov)
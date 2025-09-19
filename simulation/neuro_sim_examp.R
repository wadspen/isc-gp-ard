# Load the neuRosim package
library(neuRosim)

# 1. Define the temporal design of the fMRI experiment
# This example simulates a 200-second experiment with a TR of 2 seconds.
# It includes 20-second ON blocks starting at 1, 41, 81, 121, and 161 seconds.
# The effect size is set to 1, and a double-gamma HRF is used.
design <- simprepTemporal(totaltime = 200, 
                          onsets = seq(1, 200, 40), 
                          durations = 20, 
                          TR = 2, 
                          effectsize = 1, 
                          hrf = "double-gamma")

# 2. Define the spatial characteristics of the activated regions
# This example defines two spherical regions with specified coordinates and radii.
# The image dimensions are 64x64.
region <- simprepSpatial(regions = 2, 
                         coord = list(c(32, 15), c(57, 45)), 
                         radius = c(10, 7), 
                         form = "sphere")

# 3. Simulate the 4D fMRI data
# This simulates a 4D fMRI dataset (time x x-dim x y-dim).
# SNR is set to 1, and no noise is added for simplicity.
out <- simVOLfmri(design = design, 
                  image = region, 
                  dim = c(64, 64), 
                  SNR = 1, 
                  noise = "none")

# 4. Visualize a time series from an activated voxel
# Plot the time series of the voxel at coordinates (32, 15)
plot(out[32, 15, ], type = "l", 
     main = "BOLD Signal in Activated Voxel (32,15)", 
     xlab = "Scan Number", 
     ylab = "Signal Intensity")

# 5. Visualize a slice of the simulated fMRI data at a specific time point
# Display the 10th scan (time point) as a grayscale image.
image(1:64, 1:64, out[,,50], col = grey(0:255/255), 
      main = "Simulated fMRI Slice (Scan 10)", 
      xlab = "X-coordinate", 
      ylab = "Y-coordinate")

## Objective: To create a fully autonomous data analysis pipeline to search for Fast Radio Bursts (FRBs) working 24*7 for all observing frequency bands of Giant Metrerwave Radio Telescope, Kholad, Pune , India


![GMRT-SHM-AstroAcclerate pipeine work flow](https://github.com/Rushi314/astro-accelerate-for-GMRT-with-Shared-Memory/blob/master/NCRA%20(1).png)



## Updates to the astroacclerate pipeline for realtime processing of data from GMRT

## i) Python script documentation - spsplotii_frb_detect.py
#### Matched filtering with Gaussian pulse to identify the FRBs from the actroacclerate pipeline output

To run the code, use:

In a terminal window, type

    python spsplotii_frb_detect.py fname --ls 1.0 --auto True

In the command above, replace the “fname” with filename of the input data file. It is recommended to run the code with a low value of “--ls” such as 1.0. This will allow for maximum data capture from the input file such that the mean SNR of the data is much lower than the SNR of the FRB.

#### FRB Detection algorithm(Matched Filtering with Gaussian Pulse profile):

The input data file is read and some data is clipped from it. This clipping depends upon values of “ls”, “us”, “lt”, “ut”, “lw”, “uw”, “ld” and “ud” - these are limits upon the minimum and maximum SNR, time, pulse width, and DM respectively. Arrays which store values of DM, SNR, width, number of pulses and time are returned after clipping. The detection algorithm works on these arrays.

Since the SNR is most distinguishing feature of an FRB, the detections are made using this feature of the data. First, a grid of DM values is created (with resolution of 5). We then go to every DM grid point and apply the following formula:

![](https://latex.codecogs.com/gif.latex?f%20%3D%20e%5E%7B%5Cfrac%7B-r%5E2%7D%7B2R%5E2%7D%7D)

![](https://latex.codecogs.com/gif.latex?w%20%3D%20%5Cfrac%7B%5Csum_%7B%7D%5E%7B%7D%20s%5E3%20f%7D%7B%5Csum_%7B%7D%5E%7B%7Df%7D)

In the equations above, f is a gaussian weightage of each SNR, depending upon how far it is (r) from a Dispersion Measure grid point. R is the radius of this gaussian ball and is set at 5. SNR is represented by s in the formula above. However, we may encounter a situation where there are very many points near the FRB DM but that they may have low SNR values - being the noise. In such cases, we still want to be able to detect that an FRB signal is hidden there. Therefore, the power of 3 was chosen in the formula above such that the code is very sensitive to any pulse of a high SNR value. Simultaneously, the maximum values of width, time, and SNR are found around every grid point.

The value stored in w provides a distinction between FRBs and the background noise. Finally a threshold value is set at 1.3 and the code checks if any of w values is higher than the mean of w multiplied by this threshold. 
If there is any data that crosses this threshold, we treat it as an FRB detection. Sometimes, consecutive DM points on the grid can have overlap from FRB data. Therefore, the grid point with maximum w value is chosen among them.
Lastly, the maximum SNR value is found around this grid point and its associated DM, width and time values are printed to screen.
 

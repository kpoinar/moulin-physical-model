# Moulin Shape (MouSh) model
The model is physically based and simulates the vertical shape of a moulin over time.  MouSh represents a moulin as a stack of egg-shaped cross-sections.  The egg shape has two radii: a minor radius (R_1, smaller) and a major radius (R_2, larger).  The "egg" is a half-ellipse (radii R_1 and R_2), half-circle (radius R_1).

MouSh accepts a water input at the top of the moulin and discharges a water output through the bottom through a subglacial channel.  At each timestep, the shape of the moulin is calculated as a function of viscous deformation, elastic deformation, turbulent energy dissipation through melting, and refreezing.  At each timestep, the water level, cross-sectional area of the subglacial channel, and discharge is calculated using the popular Schoof (2010) equations.

## Main file is call_moulingeom_series_H.m
We should give this a better name.

## Troubleshooting:
### Smooth function - if "Error in importTz" appears
you need to download the function fastsmooth in matworks files exchange. It is available at https://www.mathworks.com/matlabcentral/fileexchange/19998-fast-smoothing-function. Keep it in a functions folder, making sure that the matlab path is set to include it. 

### ColorBrewer - if "Undefined function or variable 'brewermap'." appears
you need to download the function ColorBrewer in matworks files exchange. it is available at https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps. Keep it in a functions folder, making sure that the matlab path is set to include it.



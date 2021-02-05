# Moulin Shape (MouSh) model

The model simulates the shape of a moulin over time.  Moulins are often represented as straight vertical cylinders or cones, but the Moulin Shape model ([MouSh](https://www.urbandictionary.com/define.php?term=moush)) allows deviations from the cylinder.  One could think of a cylinder as a stack of circular cross-sections; MouSh represents a moulin as a stack of egg-shaped cross-sections.  The egg shape has two radii: a minor radius (R_1, smaller) and a major radius (R_2, larger).  The "egg" is a half-ellipse (radii R_1 and R_2), half-circle (radius R_1).
The radii of the "egg" vary vertically and in time.

MouSh accepts a water input at the top of the moulin and discharges a water output through the bottom through a subglacial channel.  At each timestep, the shape of the moulin is calculated as a function of viscous deformation, elastic deformation, turbulent energy dissipation through melting, and refreezing.  At each timestep, the water level, cross-sectional area of the subglacial channel, and discharge is calculated using the Schoof (2010) equations for a subglacial conduit, modified as described in Covington et al. (2020).

With parameters and melt forcings typical of the western Greenland ablation zone, the MouSh radii expand and contract on the order of 20\% daily (Andrews et al., 2020).

## Physical processes represented in MouSh


### Elastic deformation (inward/outward) of the moulin walls

See elastic.m

Elastic deformation is positive (opening) above the water line.

Elastic deformation may be positive (opening) or negative (closing) below the water line.


### Viscous deformation (inward/outward) of the moulin walls

See creep.m

Viscous deformation is positive (opening) above the water line.

Viscous deformation may be positive (opening) or negative (closing) below the water line.


### Viscous deformation of the ice column

See deformGlen.m

Shear stresses within the ice column cause ice to flow, following Glen's Flow Law.  This deforms the moulin by displacing each vertical cross-section level slightly (e.g., a shearing deck of egg-shaped cards, or a coiled, stacked Slinky being pushed).


### Melting of the submerged moulin walls via turbulent dissipation

See turbulence.m

Water velocity through the moulin (calculated according to mass conservation) and the Darcy-Weisbach equation (including variable roughness or friction factors) determines the amount of wall melt-out induced by the turbulent, net downward flow of water through the moulin.


### Melting of the moulin walls above the water line via open-channel flow dissipation

See openchannel.m

A "waterfall" sort of feature is often observed in water flowing from the ice-sheet surface into the moulin.  This function calculates the wall melt-out that results from this turbulent, open-channel flow.  


### Evolution of the subglacial channel size and water flux

See subglacialsc.m 

The moulin water level and outflow is calculated based on the Schoof (2010) equations for a subglacial channel, with modifications by Covington et al. (2020).  MouSh assumes that a single, uniform subglacial channel (length L ~ tens of km) connects the moulin to the ocean.


## Running the model 

MouSh runs in MATLAB.  A parallel version in Python will soon be available (Celia Trunz).

The main file is call_moulingeom_series_H.m
We should give this a better name.

User-defined variables include ice thickness (H), friction factors (fR), regional stresses (sigx, sigy, tauxy), ice temperature (T), shear modulus (Y), subglacial ice softness (A), and ice column creep enhancement factor (E).

Length of model run, timestep, and vertical node spacing are all variable.  Defaults are 30 days, 300 seconds, and 1 meter, respectively.

User-selected forcings include meltwater input (Qin).

Model outputs are the moulin geometry (R_1(z) and R_2(z)), water level (hw), outflux (Out), and subglacial channel size (S) at every timestep.


## Troubleshooting

### Smooth function - if "Error in importTz" appears

The function fastsmooth.m is missing. It is available at the [Mathworks File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/19998-fast-smoothing-function). Keep it in a functions folder, making sure that the MATLAB path is set to include it. 

### ColorBrewer - if "Undefined function or variable 'brewermap'." appears

The function ColorBrewer.m is missing.  It is available at the [Mathworks File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps). Keep it in a functions folder, making sure that the MATLAB path is set to include it.


# References

Andrews, Poinar, and Trunz (2020).  Physical controls on moulin geometry and evolution within the Greenland ice sheet.  Soon available on The Cryosphere Discussions.

Covington, M. D., Gulley, J. D., Trunz, C., Mejia, J., & Gadd, W. (2020). Moulin Volumes Regulate Subglacial Water Pressure on the Greenland Ice Sheet. Geophysical Research Letters, 47(20). https://doi.org/10.1029/2020GL088901

Schoof, C. (2010). Ice-sheet acceleration driven by melt supply variability. Nature, 468(7325), 803–806. https://doi.org/10.1038/nature09618



# Contact Us

Kristin Poinar, kpoinar@buffalo.edu

Lauren Andrews, lauren.c.andrews@nasa.gov

Celia Trunz, cctrunz@uark.edu


# Contributors

Lauren C. Andrews, NASA Goddard Space Flight Center

Kristin Poinar, University at Buffalo

Celia Trunz, University of Arkansas

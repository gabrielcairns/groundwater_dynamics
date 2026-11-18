Code, written in 2024, used in the article "Groundwater dynamics beneath a marine ice sheet" by Gabriel J. Cairns, Graham P. Benham and Ian J. Hewitt. The code was written in MATLAB R2024a. Data used in the codes and for the creation of plots is also included.

Many scripts used in creation of figures use colormaps from the MATLAB package:
	cmocean perceptually-uniform colormaps (Thyng et al., 2016)
In addition, the script ross_sea_basin_data.m requires the following MATLAB packages:
	Antarctic Mapping Tools (Greene et al., 2017)
	Bedmap2 Toolbox for Matlab (Fretwell et al., 2017; Greene et al., 2017)
	Antarctic boundaries, grounding line, and masks from InSAR (Rignot et al., 2013; Greene et al., 2017; Mouginot et al., 2017)
This script also uses the data sets in netcdf format:
	Ross_Embayment_basement_filt.nc
	Ross_Embayment_sediment.nc
which are available from Tankersley et al., 2022, and the dataset kingslake_gl_min.csv, which has been reverse-engineered by digitising a plot in Kingslake et al., 2018.

An index of the respective scripts and data is as follows, sorted by the order in which they are used in the article:

SCRIPTS ----------------

plot_nose_and_lens_separately.m
Creates Fig. 2.

periodic_gl.m 
Creates Fig. 3 and supplementary animations S1–3. Produces the data ICfast.mat, ICmed.mat, and ICslo.mat, which it also uses.

plot_pocket_steady_states.m
Creates Fig. 4.

pocket_dynamics.m
Creates Fig. 5 and supplementary animations S4 and S5.

pocket_gl_advance.m
Creates Fig. 6 and supplementary animation S6.

ross_sea_basin_data.m
Creates Fig. 7 and the data bHdata.mat (b and S averaged along a transect). Uses the data and packages mentioned above.

ross_steady_states.m
Creates Fig. 8. Uses the data bHdata.mat.

ross_sea_basin.m 
Creates the data ross_XXXXX.mat (solutions) and ross_XXXXX_hic.mat (initial conditions), the latter of which it also uses, where XXXXX represents one of {1e-10, 3e-12, 1e-12, 3e-13, 1e-13, 3e-14}. Uses the data bHdata.mat.

plot_ross_results.m (Figure 9-11)
Creates Figs. 9-11 and supplementary animations S7-12. Uses the data ross_XXXXX.mat, as above, and bHdata.mat.

FUNCTIONS -----

default_params.m
Initialises parameter values and scalings.

draw_frame.m
Helps with the creation of supplementary animations.

evolve_h_func.m
Objective function solved at each timestep in order to numerically evolve the freshwater-saltwater interface.

init_axes.m
Initialises axes for consistent plotting style.

SSA_ice.m
Solves the shallow-stream approximation ice sheet model used in other scripts.

viscousIce.m
Solves a general viscous ice sheet model (no longer used).

write_frame.m
Helps with the creation of supplementary animations.


DATA -------------------

/periodic_init/
Contains initial conditions h for the periodic solutions displayed in Figs. 3a-c and supplementary animations S1–3 respectively, and Supplementary Figure S1. 

bHdata.mat
Basement geometry b and S, laterally averaged along a transect, along with a curvilinear coordinate x representing distance along a transect. 

kingslake_gl_min.csv
Most retreated grounding line position, ~10 kyr before present. Obtained by digitising a plot in Kingslake et al., 2018.

/ross_results/
Periodic solutions for Ross Sea experiment. Includes data plotted in Figs. 9-11 and supplementary animations S7-12, along with Table 3.

/ross_init/
Initial conditions h for the periodic solutions above.


Bibliography ===========

Thyng, K. M., Greene, C. A., Hetland, R. D., Zimmerle, H. M., & DiMarco, S. F. (2016). True colors of oceanography: Guidelines for effective and accurate colormap selection. Oceanography, 29(3), 9-13.
Greene, C. A., Gwyther, D. E., & Blankenship, D. D. (2017). Antarctic mapping tools for MATLAB. Computers & Geosciences, 104, 151-157.
Fretwell, P., Pritchard, H. D., Vaughan, D. G., Bamber, J. L., Barrand, N. E., Bell, R., ... & Zirizzotti, A. (2013). Bedmap2: improved ice bed, surface and thickness datasets for Antarctica. The cryosphere, 7(1), 375-393.
Rignot, E., Jacobs, S., Mouginot, J., & Scheuchl, B. (2013). Ice-shelf melting around Antarctica. Science, 341(6143), 266-270.
Mouginot, J., Scheuchl, B. & Rignot, E. (2017). MEaSUREs Antarctic Boundaries for IPY 2007-2009 from Satellite Radar, Version 2 [Data Set]. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/AXE4121732AD.
Tankersley, Matthew; Horgan, Huw J; Smith Siddoway, Christine S; Caratori Tontini, Fabio; Tinto, Kirsty (2022): Basement topography and sediment thickness beneath Antarctica's Ross Ice Shelf imaged with airborne magnetic data [dataset]. PANGAEA, https://doi.org/10.1594/PANGAEA.941238
Kingslake, J., Scherer, R. P., Albrecht, T., Coenen, J., Powell, R. D., Reese, R., ... & Whitehouse, P. L. (2018). Extensive retreat and re-advance of the West Antarctic Ice Sheet during the Holocene. Nature, 558(7710), 430-434.



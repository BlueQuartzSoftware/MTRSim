MTRsim
-----------------
Copyright (c) 2025 US Govt

Welcome to the MTRsim software.  This software contains code in the form
of *.m files to generate synthetic volumes of microstructures containing
microtexture regions. The inputs are a collection of sub-/component-ODFs
and their corresponding volume fractions and spatial parameters, and the
dimensions of the output voxel-based rectangular prism.  The output is a
simulated random draw of a volume containing random regions of similar 
crystallographic orientation. 

Distribution A. Approved For Public Release: Distribution Unlimited. AFRL-2025-1631

Disclaimer
-----------------
The software is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE 
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

Feedback
-----------------
Daniel.Sparkman.1 at us dot af dot mil


USAGE
-----------------
Run the simulation from the file: 
simulate_MTRs.m

INPUTS
-----------------
###The simulation volume parameters are set in simulate_MTRs.m
xlen   : length of volume in x
ylen   : length of volume in y
zlen   : length of volume in z (can be set to zero to create a 2D simulation)
dx     : voxel spacing in x
dy     : voxel spacing in y
dz     : voxel spacing in z

### The MTR spatial parameters are set in simulate_MTRs.m

volume_fractions   :  volume fraction of each MTR "class"
theta_list         :  correlation lengths of plurigaussian model
nugvar             :  nugget variance (not used)

### The plurigaussian model is used to simulate the underlying MTR assignment for each voxel 

### The plurigaussian parameters are set in PGRF_simulation.m

correlation_function_selected   : generally set to 'anisotropic'
corr_func_name                  : correlation model, set to 'exp' or 'gaussian'
boundary_conditions             : 'periodic' or 'nonperiodic'
mean_function_selected          : generally set to 'stationary'

Based on the MTR assignments, for each voxel a random crystallographic 
orientation is drawn from the corresponding MTR-/component-ODF.
The finite mixture model for the ODF is essentially a collection of 
normalized histograms over the Euler angles, contained in simulation_ODF.mat

- Parameterized families of models are not provided
ODF_best                        : typically components are estimated from data 


File contents and structure
-----------------
simulate_MTRs.m
PGRF_simulation.m
select_AR.m
norminv_code.m
qsimvn.m
gp_generator.m
calc_ODF.m
convert_ODF_to_PF.m
eval_AR.m
sample_N_orientations_from_ODF.m
view_IPF_map.m
unit_triangle_IPF_coords.m
generate_inverse_pole_figure.m
sample_orientation_from_ODF.m
symmetric_euler_angles.m
IPF_colors.m
plot_AR.m
plot_PF.m
load_MTR_data.m
blank_ODF.mat
PF_coords_for_ODF.mat
PF_indexed.mat
simulation_ODF.mat
uniformODF.mat
uniformPF.mat
README.txt

Acknowledgements
-----------------
The qsimvn.m provided is based on the function of the same name by Alan Genz. It has been
modified slightly for use with this software (MTRsim).  In the file qsimvn.m
the full acknowledgement as well as original redistribution conditions, 
copyright notice, and disclaimer are contained. In accordance with the redistribution
conditions, no claim of endorsement by Alan Genz is made for MTRsim.  

Citing
-----------------
The software has not been published yet. Future citations would be formatted 
according to the following:
The appropriate citations are:
@article{MTRsim, 
 title={Spatial Statistical Simulation of Microtexture},
 author={Daniel Sparkman, Laura Homa, Adam Pilchak, and Harry Millwater},
 journal={###},
 volume={#},
 number={#},
 month={###.},
 year={2024}
}
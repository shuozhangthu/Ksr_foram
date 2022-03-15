Model and data analysis scripts for "Consistent Sr partitioning in foraminiferal and inorganic calcite".
 
This directory includes:
-Model input data.
   (1) Core-top data from Yu et al. (2014) for C.wuellerstorfi and C.mundulus. Files: carbonate_parameters.mat, carbonate_parameters_Cm.mat
   (2) Script for calculating kinetic (attachment/detachment frequencies) and thermodynamic parameters. File: kinetic_thermodynamic_parameters.m

-Box model scripts for C.wuellerstorfi. 
   (1) Scripts for calculating stable calcifying fluid  composition at varying seawater depths in four oceans. Files: function_chen_depth_a.m, function_chen_depth_i.m, function_chen_depth_n.m, function_chen_depth_p.m
   (2) Calculation of model output. Includes the solution chemistry of seawater in four Oceans and the C.wuellerstorfi calcifying fluid, modeled calcite precipitation rate and the foraminiferal KSr. Files: run_chen_depth_a.m, run_chen_depth_i.m, run_chen_depth_p.m, run_chen_depth_n.m
   (3) Plot the figures of Fig.3, Fig.4 and Fig.S1. Files: figure_sw_cf_Rp_depth.m, figure_sw_cf_Rp_depth2.m, figure_carbonate_depth.m, figure_ksr_omega.m

-Box model scripts for C.mundulus.
   (1) Script for calculating stable calcifying fluid  composition at varying seawater depths in the Atlantic ocean. File: function_Cm_A.m
   (2) Calculation of model output. Includes the solution chemistry of seawater in the Atlantic Ocean and the C.mundulus calcifying fluid, modeled calcite precipitation rate and the foraminiferal KSr. File: run_Cm_A.m      

Note: File [equic.m] are used to calculate seawater carbonate chemistry in function files. File [solveP.m] are used to calculate the probability of each kink site in run files and function files.

Reference:
Yu, J., Elderfield, H., Jin, Z., Tomascak, P. and Rohling, E.J. (2014) Controls on Sr/Ca in benthic foraminifera and implications for seawater Sr/Ca during the late Pleistocene. Quaternary Science Reviews 98, 1-6.
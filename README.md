# mo_lons

Matlab codebase for multi-objective LONs. 

The subdirectory gecco_2019 contains the codebase for the paper

Jonathan E. Fieldsend and Khulood Alyahya. 2019. 
Visualising the Landscape of Multi-Objective Problems using Local Optima Networks. 
In Genetic and Evolutionary Computation Conference Companion (GECCO ’19 Companion), 
July 13–17, 2019, Prague, Czech Republic. ACM, New York, NY, USA,

Paper available at

Institutional repository: http://hdl.handle.net/10871/36891

Publisher DOI: https://doi.org/10.1145/3319619.3326838

Note the gecco_2019 codebase is illustrative of proceedures and not optimised for speed!

FUNCTION OVERVIEW

generate_GECCO_2019_plots.m --- script will regenerate all the PLON, DNON and PLOS-net plots from the paper

gecco_workshop_2019_problem1.m, gecco_workshop_2019_problem2.m and gecco_workshop_2019_problem3.m --- cost functions used in illustrations in the paper. The last takes an argument holding the grid for each objective (this data is stored in the two .mat files loaded in regenerate_GECCO_2019_plots.m).

exhaustive_generate_lon -- searches over the (discretised) planar cost function exhaustively at the given resolution, and restuns matrices and structures used by the three network generators

plo_PLOS_net.m -- plots the PLOS-net given the arguments

process_d_lon.m -- generates the DNON matrices

process_p_lon.m -- generates the PLON matrices



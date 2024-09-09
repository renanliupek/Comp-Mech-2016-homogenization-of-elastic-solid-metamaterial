# Comp-Mech-2016-homogenization-of-elastic-solid-metamaterial

matlab script implementing the computational homogenization of locally resonant solid metamaterials with an elastic behavior, as in the article below.

[1] Sridhar, Ashwin, Varvara G. Kouznetsova, and Marc GD Geers. "Homogenization of locally resonant acoustic metamaterials towards an emergent enriched continuum." Computational mechanics 57 (2016): 423-435.

## Note

1- run `start.m` to include the the path to the folders *fun*, *tensorlab* and *meshes*.

2- export comsol 5.4 mesh file `.mphtxt` of your favorite geometry to the folder *meshes*, or use the mesh file provided.

3- FATAL: remember to rename the exported comsol mesh `.mphtxt` to `.txt`.


## `main_pbc.m` file description

1- `main_pbc.m` computes the homogenised material coefficients using 'periodic fluctuation' boundary conditions as in the paper.


## `main_kbc.m` file description

1- `main_kbc.m` computes the homogenised material coefficients using 'zero fluctuation' boundary conditions. This is not done in the paper but it is left as an option to the user given the popularity of this boundary condition in the computational homogenization community.

## *fun* folder
Developed functions to assist the computations in main files.

## *tensorlab* folder
Collection of functions and classes that define a tensor object.

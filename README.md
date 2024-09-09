MATLAB scripts of the paper:

Sridhar, Ashwin, Varvara G. Kouznetsova, and Marc GD Geers. "Homogenization of locally resonant acoustic metamaterials towards an emergent enriched continuum." Computational mechanics 57 (2016): 423-435.

## Note

1- run `start.m` to include the the path to the folders *fun*, *tensorlab* and *meshes*.

2- export comsol 5.4 mesh file `.mphtxt` of your favorite geometry to the folder *meshes*, or use the mesh file provided.

3- FATAL: remember to rename the exported comsol mesh `.mphtxt` to `.txt`.

## `main_pbc.m` file description

1- `main_pbc.m` computes the homogenized material coefficients using 'periodic fluctuation' boundary conditions for the unit cell as in the paper.

2-  The input for `main.m` are: the mesh file in the folder *meshes* with the list of tags of the corresponding comsol model, and the material properties of the elastic solid phases which can be defined within the script `main.m`.

3-  The appropriate mesh is created by exporting a text file `.mphtxt` from comsol 5.4, see 'Note' 2nd and 3rd item.

## `main_kbc.m` file description

1- `main_kbc.m` computes the homogenized material coefficients using 'zero fluctuation' boundary conditions for the unit cell as in the paper. This is a popular boundary condition in computational homogenization community. If the representative volume element (defined as the unit cell) is too small, then it does not approximate well the overall material properties.

2-  The input for `main.m` are: the mesh file in the folder *meshes* with the list of tags of the corresponding comsol model, and the material properties of the elastic solid phases which can be defined within the script `main.m`.

3-  The appropriate mesh is created by exporting a text file `.mphtxt` from comsol 5.4, see 'Note' 2nd and 3rd item.

## *fun* folder
Developed functions to assist the computations in `main.m`.

## *tensorlab* folder
Collection of functions and classes that define a tensor object.

## *meshes* folder
It contains the mesh files of the unit cell shown in the article or any other unit cell you may want to homogenise. 


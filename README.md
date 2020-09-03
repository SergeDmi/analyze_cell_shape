# Analyze Cell Shape
This is a software suite to analyze the shape of triangulated objects stored, into ply files.  
This is mostly designed for rod-like cells, e.g. yeast, but any closed triangulated mesh stored in a ply file should work.  
To segment shapes from microscopy images, we recommend Limeseg : https://imagej.net/LimeSeg

## Operation
We provide a set of simple functions to analyze the shape of a single cell, or to compare cells at different stages. The simplest is to go through the interface of *analyze_shape(shape)* a function that takes as argument the shape of a cell (see below : I/O variables). This outputs *results*, a structure containing a set of results. To compare several stages for a set of cells, the function *summary_experiments* should be used. It takes for input *experiments*, a formated list of severall cell shape shapes for different stages (see below : I/O variables). Several simple examples of usage are provided.

### Single cell analysis.
See demo_single_cell.m  

### Comparing states
To compare the states of several cells of a number for a various number of stages (i.e. each cell has different states), see :  
- demo_compare_2_stages.m  
- demo_compare_4_stages.m  

## Input/output variables
### Shape
Shape is the input argument for *analyze_shape*.  A shape is a structure containing *points* - an array of vertices -, *faces* - the array of triangles -, *normals* - an array of the normals to the faces-.  
- shape.points : $ 3 \times N_v$ array with $N_v$ the number of vertices. Each element is a float representing a position.  
- shape.faces : $ 3  \times N_f$ array with $N_f$ the number of faces. Each element is an float representing a vertex index.  
- shape.normals : $ 3 \times N_f$ array. Each element is a float representing an element of the normal vector.  
A *shape* is the output from *import_ply*, itself a wrapper for *ply_read.m*, and external module.

### Experiment
An *experiment* is a structure containing several states. Each state is a structure containing a *name* (usually the filename) and a *shape*. E.g. :   
- experiment.states(1).name='cell_1_stage_1.ply'  
- experiment.states(1).shape=import_ply('cell_1_stage_1.ply')  
- experiment.states(2).name='cell_1_stage_2.ply'  
- experiment.states(2).shape=import_ply('cell_1_stage_2.ply')  

Experiment can also store an *analysis*, the results of *analyze_shape* :  
- experiment.states(1).analysis=analyze_shape(experiment.states(1).shape)  

### Experiments
Experiments is the input argument for *summary_experiments*. Experiment is a list of several experiments :  
- experiments(1).states(1).name='cell_1_stage_1.ply'  
- experiments(1).states(1).shape=import_ply('cell_1_stage_1.ply')  
- experiments(2).states(1).name='cell_2_stage_1.ply'  
- experiments(2).states(1).shape=import_ply('cell_2_stage_1.ply')  

Each element of experiments may also store an *analysis*, the result of *analyze_shape* :  
- experiments(i).states(j).analysis=analyze_shape(experiments(i).states(j).shape)  

### Analysis
An *analysis* is the results of *analyze_shape* analyzing a *shape*. It contains many derived quantities. See the code for the definition of these quantities. The most important quantities are :   
- *analysis.length* : the length of the shape  
- *analysis.volume* : the volume of the shape  
- *analysis.surface* : the surface of the shape  
- *analysis.mean_circ* : mean circularity of slices in the cell, weighted by surface area  
- *analysis.mean_pers* : mean perimeter of slices in the cell, weighted by surface area  
- *analysis.mean_area* : mean area of slices in the cell, weighted by surface area  


# Serge Dmitrieff -- http://biophysics.fr

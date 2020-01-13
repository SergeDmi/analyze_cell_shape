% Loading local file directory
addpath(genpath('.'));
% Adding names
cell=import_ply('sphere.ply');
analysis=analyze_shape(cell)

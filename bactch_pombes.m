%% We're gonna do a batch analysis !
% We scan over many folders (experiments)
% for which we get several sub-experiment
% for which there is pre and post
% for which there are several cells
% ...
% Let's go !

%% Ok now we have to declare the root and folders
root='pombe_3D_files/';
folders={'180214/','180220/','180222/'};
[experiments]=load_exps(root,folders);

[res,experiments]=summary_experiments(experiments);
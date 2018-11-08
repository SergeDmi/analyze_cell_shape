function [ experiment ] = align_experiments( experiment,options )
% We align the two cells
%   Nothing magical

%% We  first test the input
if nargin<2
    options=pombe_default_options();
end
if isfield(options,'pixel_size')
	pixel_size=options.pixel_size;
else
	disp('Warning : could not find pixel_size in options');
	pixel_size=1;
end
%% Loading and centering points
pp1=experiment.pre_pombe.points*pixel_size;
pp2=experiment.post_pombe.points*pixel_size;
c1=mean(pp1);
c2=mean(pp2);
pp1=pp1-ones(size(pp1,1),1)*c1;
pp2=pp2-ones(size(pp2,1),1)*c2;

%% PCA to find the main axis of pre points
[co1,p1]=princom(pp1);
p2=pp2*co1;

if options.aligning > 0
	 [~,p2,~]=princom(pp2);
end

%% Getting away with murder
experiment.pre_pombe.points=p1/pixel_size;
experiment.post_pombe.points=p2/pixel_size;


[experiment.pre_pombe.backbone,experiment.pre_pombe.fitted_backbone,experiment.pre_pombe.fitted_arclen]=get_rolling_backbone(p1,options.backbone);
[experiment.post_pombe.backbone,experiment.post_pombe.fitted_backbone, experiment.post_pombe.fitted_arclen]=get_rolling_backbone(p2,options.backbone);

end

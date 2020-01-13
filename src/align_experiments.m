function [ experiment ] = align_experiments( experiment,options )
% We align the two cells
%   Nothing magical

%% We  first test the input
if nargin<2
    analysis_options=analysis_default_options();
    summary_options=summary_default_options();
else
  analysis_options=options.analysis_options;
  summary_options=options.summary_options;
end
if isfield(options,'pixel_size')
  	pixel_size=options.pixel_size;
else
  	disp('Warning : could not find pixel_size in options');
  	pixel_size=1;
end
% Number of different states/stages
n_states=numel(experiment.states)


%% Loading and centering points
% First, the first cell
pp1=experiment.states(1).shape.points*pixel_size;
% centering
c1=mean(pp1);
pp1=pp1-ones(size(pp1,1),1)*c1;
% aligning
[co1,p1]=get_pca_cov(pp1,1);
experiment.states(1).shape.points=p1/pixel_size;

for i=2:n_states
  pp2=experiment.states(i).shape.points*pixel_size;
  c2=mean(pp2);
  pp2=pp2-ones(size(pp2,1),1)*c2;
  if summary_options.align_to_first
    % We align to the first stage
    p2=pp2*co1;
  else
    [~,p2,~]=get_pca_cov(pp2,1);
  end
  experiment.states(i).shape.points=p2/pixel_size;

end

% Legacy code
% backbone computation should be in analysis, not alignment !
%[experiment.pre_pombe.backbone,experiment.pre_pombe.fitted_backbone,experiment.pre_pombe.fitted_arclen]=get_rolling_backbone(p1,options.backbone);
%[experiment.post_pombe.backbone,experiment.post_pombe.fitted_backbone, experiment.post_pombe.fitted_arclen]=get_rolling_backbone(p2,options.backbone);

end

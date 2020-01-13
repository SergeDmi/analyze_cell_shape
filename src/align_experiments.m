function [ experiment ] = align_experiments( experiment,options )
% We align the two cells
%   Nothing magical

%% We  first test the input
if nargin<2
    options=summary_default_options();
end
% Number of different states/stages
n_states=numel(experiment.states);


%% Loading and centering points
% First, the first cell
pp1=experiment.states(1).shape.points;
% centering
c1=mean(pp1);
pp1=pp1-ones(size(pp1,1),1)*c1;
% aligning
%[co1,p1]=get_pca_cov(pp1,1);
[co1,p1]=princom(pp1);

for i=1:n_states
  pp2=experiment.states(i).shape.points;
  c2=mean(pp2);
  pp2=pp2-ones(size(pp2,1),1)*c2;
  if options.align_to_first
    % We align to the first stage
    p2=pp2*co1;
  else
    %[~,p2,~]=get_pca_cov(pp2,1);
    [~,p2,~]=princom(pp2);
  end
  experiment.states(i).shape.points=p2;

end

% Legacy code
% backbone computation should be in analysis, not alignment !
%[experiment.pre_pombe.backbone,experiment.pre_pombe.fitted_backbone,experiment.pre_pombe.fitted_arclen]=get_rolling_backbone(p1,options.backbone);
%[experiment.post_pombe.backbone,experiment.post_pombe.fitted_backbone, experiment.post_pombe.fitted_arclen]=get_rolling_backbone(p2,options.backbone);

end

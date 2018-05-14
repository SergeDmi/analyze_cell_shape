function [ experiment ] = align_experiments( experiment )
% We align the two cells 
%   Nothing magical

%% Loading and centering points
pp1=experiment.pre_pombe.points;
pp2=experiment.post_pombe.points;
c1=mean(pp1);
c2=mean(pp2);
pp1=pp1-ones(size(pp1,1),1)*c1;
pp2=pp2-ones(size(pp2,1),1)*c2;

%% PCA to find the main axis of pre points
[co1,p1]=princom(pp1);
p2=pp2*co1;

%% Getting away with murder
experiment.pre_pombe.points=p1;
experiment.post_pombe.points=p2;

end


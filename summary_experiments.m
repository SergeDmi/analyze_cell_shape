function [ results ] = summary_experiments( experiments )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nexp=numel(experiments);
aspect=zeros(Nexp,2);
surface=zeros(Nexp,2);
volume=zeros(Nexp,2);
length=zeros(Nexp,2);
radius=zeros(Nexp,2);
interest=true(Nexp,1);
for n=1:Nexp
	surface(n,1)=experiments(n).pre_analysis.surface;
	volume(n,1)=experiments(n).pre_analysis.volume;
	aspect(n,1)=experiments(n).pre_analysis.central_circularity;
	length(n,1)=experiments(n).pre_analysis.length;
	
	surface(n,2)=experiments(n).post_analysis.surface;
	volume(n,2)=experiments(n).post_analysis.volume;
	aspect(n,2)=experiments(n).post_analysis.central_circularity;
	length(n,2)=experiments(n).post_analysis.length;	
end

results.surface=surface;
results.volume=volume;
results.aspect=aspect;
results.length=length;

results.volume_change=(diff(volume,1,2));
results.surface_change=(diff(surface,1,2));
results.aspect_change=(diff(aspect,1,2));
results.length_change=(diff(length,1,2));

[~,nv]=min(results.volume_change);
interest(nv)=true;
experiments(nv).tag='Minimal volume difference';

[~,nv]=max(results.volume_change);
interest(nv)=true;
experiments(nv).tag='Maximal volume difference';

[~,nv]=max(results.length_change);
interest(nv)=true;
experiments(nv).tag='Maximal length difference';

[~,nv]=min(results.aspect_change);
interest(nv)=true;
experiments(nv).tag='Maximal aspect change';
%results.radius_change=(diff(radius,1,2));

interesting=experiments(interest);
N=numel(interesting);
for n=1:N
	[~,points_pre] =princom(interesting(n).pre_pombe.points);
	[~,points_post]=princom(interesting(n).post_pombe.points);
	points_post(:,2)=points_post(:,2)+100;
	figure
    hold all
	if isfield(interesting(n),'tag')
		if ~isempty(interesting(n).tag)
			title([interesting(n).tag ' for n=' num2str(n)]);
		else
			title(num2str(n));
		end
	else
		title(num2str(n));
	end
    scatter3(points_pre(:,1),points_pre(:,2),points_pre(:,3),5,'k')
	scatter3(points_post(:,1),points_post(:,2),points_post(:,3),5,'r')
	axis equal;
end
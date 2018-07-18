function [ results,experiments ] = summary_experiments( experiments,options )
% runs over experiments
%   Detailed explanation goes here


if nargin<2
    options=pombe_default_options();
end

Nexp=numel(experiments);
aspect=zeros(Nexp,2);
surface=zeros(Nexp,2);
volume=zeros(Nexp,2);
length=zeros(Nexp,2);
radius=zeros(Nexp,2);
interest=true(Nexp,1);

%res(Nexp).pre_per=[0,0];
pre_pers=[];
change_pers=[];
pre_ares=[];
change_ares=[];
pre_circ=[];
change_circ=[];
dist_circ=[];
for n=1:Nexp
	experiments(n)=align_experiments(experiments(n));

	if ~options.check_pairs
		checked=1;
	else
		checked=check_pair(experiments(n));
	end

	if checked>0
		%experiments(n)=compare_pombes(experiments(n));
		experiments(n).pre_analysis=analyze_pombe(experiments(n).pre_pombe);
		experiments(n).post_analysis=analyze_pombe(experiments(n).post_pombe);

		experiments(n).pre_analysis=analyze_pombe(experiments(n).pre_pombe);
		experiments(n).post_analysis=analyze_pombe(experiments(n).post_pombe);


		surface(n,1)=experiments(n).pre_analysis.surface;
		volume(n,1)=experiments(n).pre_analysis.volume;
		aspect(n,1)=experiments(n).pre_analysis.central_circularity;
		length(n,1)=experiments(n).pre_analysis.length;

		surface(n,2)=experiments(n).post_analysis.surface;
		volume(n,2)=experiments(n).post_analysis.volume;
		aspect(n,2)=experiments(n).post_analysis.central_circularity;
		length(n,2)=experiments(n).post_analysis.length;

    if (options.do_slice_analysis )
      res=compare_pombes(experiments(n));
  		pre_pers=[pre_pers res.pre_per];
  		change_pers=[change_pers res.post_per-res.pre_per];
  		pre_ares=[pre_ares res.pre_area];
  		change_ares=[change_ares res.post_area-res.pre_area];

  		pre_circ=[pre_circ res.pre_circ];
  		change_circ=[change_circ res.post_circ-res.pre_circ];
  		dist_circ=[dist_circ res.dist_circ];
    end
	end
end

if (options.do_slice_analysis )
  results.pre_pers=pre_pers;
  results.change_pers=change_pers;

  results.pre_ares=pre_ares;
  results.change_ares=change_ares;

  results.pre_circ=pre_circ;
  results.change_circ=change_circ;
  results.dist_circ=dist_circ;


  if 1
     figure
     scatter(pre_pers,pre_pers+change_pers)
     figure
     scatter(pre_ares,pre_ares+change_ares)
  end

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


if 0
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

end
end

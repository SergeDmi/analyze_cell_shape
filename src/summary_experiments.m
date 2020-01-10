function [ results,experiments ] = summary_experiments( experiments,options)
% Compare diffent states for a set of experiments
%   e.g. : several cells are considered as several experiments
%   e.g. : several time points of the same cell are several states
%      experiments is a container for several states
%      each state is a container for points, filename, analysis
if nargin<2
    analysis_options=analysis_default_options();
    summary_options=summary_default_options();
else
  analysis_options=options.analysis_options;
  summary_options=options.summary_options;
end

% Bits of counting
Nexp=numel(experiments);
n_states=numel(experiments(1).states);

% What to compare
comparisons=summary_options.comparisons;
n_comp=numel(comparisons);
results(n_states).done=0;

for s=1:n_states
  for c=1:n_comp
    results(s).(comparisons{c}).values=zeros(Nexp,1);
  end
end


%aspect=zeros(Nexp,2);
%circ=zeros(Nexp,2);
%mean_circ=zeros(Nexp,2);
%mean_area=zeros(Nexp,2);
%mean_pers=zeros(Nexp,2);
%surface=zeros(Nexp,2);
%shape_perimeter=zeros(Nexp,2);
%shape_backbone=zeros(Nexp,2);
%volume=zeros(Nexp,2);
%length=zeros(Nexp,2);
%radius=zeros(Nexp,2);
%interest=true(Nexp,1);

%res(Nexp).pre_per=[0,0];
%pre_pers=[];
%change_pers=[];
%pre_ares=[];
%change_ares=[];
%pre_circ=[];
%change_circ=[];
%dist_circ=[];


%% Looping over all experiments
for n=1:Nexp
  %% We align the shapes
	experiments(n)=align_experiments(experiments(n),options);
  checked=1;

  if summary_options.check_states
    disp('Manually checking experiments')
    disp(' press [SPACE] to keep, [q] to discard')
    checked=check_states(experiments(n));
  end
  %% Legacy code : pair checking
	%if ~options.check_pairs
	%	checked=1;
	%else
  %  disp('Manually checking experiments')
  %  disp(' press [SPACE] to keep, [q] to discard')
	%	checked=check_pair(experiments(n));
	% end

	if checked>0
		%experiments(n)=compare_shapes(experiments(n));
    for s=1:n_states
      % For each state, we analyze the shape
      experiments(n).states(s).analysis=analyze_shape(experiments(n).states(s).shape,analysis_options);

      for c=1:n_comp
        % For each state, we store what we are supposed to compare between stages
        results(s).(comparisons{c}).values(n)=experiments(n).states(s).analysis.(comparisons{c});
      end
    end
  end
	%experiments(n).pre_analysis=analyze_shape(experiments(n).pre_shape);
	%experiments(n).post_analysis=analyze_shape(experiments(n).post_states);


	%surface(n,1)=experiments(n).pre_analysis.surface;
	%volume(n,1)=experiments(n).pre_analysis.volume;
  %area(n,1)=experiments(n).pre_analysis.central_area;
  %perimeter(n,1)=experiments(n).pre_analysis.central_perimeter;
	%aspect(n,1)=experiments(n).pre_analysis.central_circularity;
	%circ(n,1)=experiments(n).pre_analysis.circularity_fiji;
	%length(n,1)=experiments(n).pre_analysis.length;
	%shape_perimeter(n,1)=experiments(n).pre_analysis.perimeter_angular_deviation;
	%shape_backbone(n,1)=experiments(n).pre_analysis.backbone_curv(1);

	%surface(n,2)=experiments(n).post_analysis.surface;
	%volume(n,2)=experiments(n).post_analysis.volume;
  %area(n,2)=experiments(n).post_analysis.central_area;
  %perimeter(n,2)=experiments(n).post_analysis.central_perimeter;
	%aspect(n,2)=experiments(n).post_analysis.central_circularity;
	%circ(n,2)=experiments(n).post_analysis.circularity_fiji;
	%length(n,2)=experiments(n).post_analysis.length;
	%shape_perimeter(n,2)=experiments(n).post_analysis.perimeter_angular_deviation;
	%shape_backbone(n,2)=experiments(n).post_analysis.backbone_curv(1);
	%experiments(n).post_analysis.backbone_curv(1)

  %% Lecacy code kept in case
  %if 0
  %  if (analysis_options.do_slice_analysis )
  %		res=compare_statess(experiments(n));
  %  		pre_pers=[pre_pers res.pre_per];
  %  		change_pers=[change_pers res.post_per-res.pre_per];
  %  		pre_ares=[pre_ares res.pre_area];
  %  		change_ares=[change_ares res.post_area-res.pre_area];

  %  		pre_circ=[pre_circ res.pre_circ];
  %  		change_circ=[change_circ res.post_circ-res.pre_circ];
  %  		dist_circ=[dist_circ res.dist_circ];

  %      mean_circ(n,1)=experiments(n).pre_analysis.mean_circ;
  %      mean_circ(n,2)=experiments(n).post_analysis.mean_circ;
  %	     mean_area(n,1)=experiments(n).pre_analysis.mean_area;
  %      mean_area(n,2)=experiments(n).post_analysis.mean_area;
  %	     mean_pers(n,1)=experiments(n).pre_analysis.mean_pers;
  %      mean_pers(n,2)=experiments(n).post_analysis.mean_pers;
  %  end
  %end
  % // end legacy code //


end



%% Wrap up :
for s=1:n_states
  for c=1:n_comp
    % For each state, we compute mean & std of
    % what we are supposed to compare between stages
    results(s).(comparisons{c}).mean=mean(results(s).(comparisons{c}).values);
    results(s).(comparisons{c}).std = std(results(s).(comparisons{c}).values);
  end
end

%% Lecacy code kept just in case
%if 0
%  if (analysis_options.do_slice_analysis )
%    results.pre_pers=pre_pers;
%    results.change_pers=change_pers;
%
%    results.pre_ares=pre_ares;
%    results.change_ares=change_ares;

%    results.pre_circ=pre_circ;
%    results.change_circ=change_circ;
%    results.dist_circ=dist_circ;


%    if summary_options.verbose>0
%       figure
%       scatter(pre_pers,pre_pers+change_pers)
%       xlabel('Intial perimenter')
%       ylabel('Final perimeter')
%       figure
%       scatter(pre_ares,pre_ares+change_ares)
%       xlabel('Intial area')
%       ylabel('Final area')
%    end
%  end
%end
% // End legacy code //



%results.surface=surface;
%results.volume=volume;
%results.perimeter=perimeter;
%results.area=area;
%results.aspect=aspect;
%results.circ=circ;
%results.length=length;
%results.shape_perimeter=shape_perimeter;
%results.shape_backbone=shape_backbone;
%results.mean_circ=mean_circ;
%results.mean_pers=mean_pers;
%results.mean_area=mean_area;

%results.volume_change=(diff(volume,1,2));
%results.surface_change=(diff(surface,1,2));
%results.aspect_change=(diff(circ,1,2));
%results.circ_change=(diff(aspect,1,2));
%results.length_change=(diff(length,1,2));
%results.backbone_recovery=(shape_backbone(:,1)-shape_backbone(:,2))./shape_backbone(:,1);
%results.perimeter_recovery=(shape_perimeter(:,1)-shape_perimeter(:,2))./shape_perimeter(:,1);

%% Legacy code
if 0
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


  if verbose>0
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
% // end legacy code //


end

function [ states] = compare_pombes( experiment,options)
% graphically compares the shape of two cells
%  experiment is a container for cell points and other information
%
% Serge Dmitrieff, IJM 2018
% www.biophysics.fr

%% Deprecated :
% Experiment now has a different structure
% Anyway, this was not put to much use


if nargin < 1
    error('First argument should be an experiment');
end

if nargin < 2
    options=analysis_default_options();
end

if ~isfield(options,'colors')
  options.colors={'k','b','r','g','m','c'}
end

pixel_size=1.0;
if isfield(options,'pixel_size')
  pixel_size=options.pixel_size;
end

%% Here we plot
n_states=numel(experiment.states);


%coeffs=experiment.pre_analysis.coeffs;
states(n_states).min=0;
mins_pos=zeros(n_states,1);
maxs_pos=zeros(n_states,1);

mins_pers=zeros(n_states,1);
maxs_pers=zeros(n_states,1);
mins_ares=zeros(n_states,1);
maxs_ares=zeros(n_states,1);

for s=1:n_states
    states(s).p1=experiment.states(s).shape.points;
    states(s).res=experiment.states(s).analysis.slices;
    mins_pos(s)=min(states(s).res.pos);
    maxs_pos(s)=max(states(s).res.pos);
    mins_pers(s)=min(states(s).res.pers);
    maxs_pers(s)=max(states(s).res.pers);
    mins_ares(s)=min(states(s).res.ares);
    maxs_ares(s)=max(states(s).res.ares);
end

hh=(min(mins_pos):options.slice.spline_dh:max(maxs_pos));

figure

for s=1:n_states
  res=states(s).res;
  p1=states(s).p1;
  res.pesp=spline(res.pos,res.pers,hh);

  res.area=spline(res.pos,res.ares,hh);

  res.circ=spline(res.pos,res.circ,hh);

  res.dist_circ=spline(res.pos,res.dist,hh);
  states(s).res=res;
  if options.verbose>0

    color=get_color(s,options.colors)
    subplot(3,n_states,s)
    ylabel('perimeter')
    hold all
    scatter(res.pos,res.pers,color)
    plot(hh,res.pesp,color)
    axis([min(mins_pos) max(maxs_pos) 0 max(maxs_pers)]);


    subplot(3,n_states,n_states+s)
    ylabel('area')
    hold all
    scatter(res.pos,res.ares,color)
    plot(hh,res.area,color)
    scatter(res.pos,res.ares,color)
     plot(hh,res.area,color)
    %axis([min([preres.pos;post.res.pos]) max([pre.res.pos;post.res.pos]) 0 max([pre.res.ares;post.res.ares])]);
    axis([min(mins_pos) max(maxs_pos) 0 max(maxs_ares)]);

    %title(titles);
    subplot(3,n_states,[(2*n_states+1):3*n_states])
    hold all
    scatter3(p1(:,1),p1(:,2),p1(:,3),5,color)
    for i=1:res.nslice
      %res.per
      %size(res.per(i).points)
        plot3(res.per(i).points(:,1),res.per(i).points(:,2),res.per(i).points(:,3),color,'LineWidth',2)
    end

    axis auto;
    axis equal;
    curr_name=experiment.states(1).name;
    if length(curr_name)>20
      curr_name=['...' curr_name(end-20:end)]
    end
    title(curr_name, 'Interpreter', 'none');
  end

end



end

function color=get_color(n,colors)
  nc=numel(colors);
  while n>nc
    n=n-nc;
  end
  color=colors{n};
end

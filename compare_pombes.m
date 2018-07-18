function [ res] = compare_pombes( experiment,options)
% graphically compares the shape of two cells
%  experiment is a container for cell points and other information
%
% Serge Dmitrieff, IJM 2018
% www.biophysics.fr

if nargin<2
    options=pombe_default_options();
end

%coeffs=experiment.pre_analysis.coeffs;


p1=experiment.pre_pombe.points;
p2=experiment.post_pombe.points;

pre.res=experiment.pre_analysis.slices;
post.res=experiment.post_analysis.slices;

hh=(min([pre.res.pos;post.res.pos])):options.spline_dh:(max([pre.res.pos;post.res.pos]));
res.pre_per=spline(pre.res.pos,pre.res.pers,hh);
res.post_per=spline(post.res.pos,post.res.pers,hh);

res.pre_area=spline(pre.res.pos,pre.res.ares,hh);
res.post_area=spline(post.res.pos,post.res.ares,hh);

res.pre_circ=spline(pre.res.pos,pre.res.circ,hh);
res.post_circ=spline(post.res.pos,post.res.circ,hh);
res.dist_circ=spline(pre.res.pos,pre.res.dist,hh);

if options.verbose>0
   figure
   subplot(3,2,1)
    hold all
    scatter(pre.res.pos,pre.res.pers,'k')
    plot(hh,res.pre_per,'k')
    scatter(post.res.pos,post.res.pers,'r')
    plot(hh,res.post_per,'r')
    axis([min([pre.res.pos;post.res.pos]) max([pre.res.pos;post.res.pos]) 0 max([pre.res.pers;post.res.pers])]);
    ylabel('perimeter')
    subplot(3,2,2)
    hold all
    scatter(pre.res.pos,pre.res.ares,'k')
    plot(hh,res.pre_area,'k')
    scatter(post.res.pos,post.res.ares,'r')
     plot(hh,res.post_area,'r')
    axis([min([pre.res.pos;post.res.pos]) max([pre.res.pos;post.res.pos]) 0 max([pre.res.ares;post.res.ares])]);
    ylabel('area')
    %title(titles);
    subplot(3,2,[3,6])
    hold all
    scatter3(p1(:,1),p1(:,2),p1(:,3),5,'k')
    for i=1:pre.res.nslice
        plot3(pre.res.per(i).points(:,1),pre.res.per(i).points(:,2),pre.res.per(i).points(:,3),'k','LineWidth',2)
    end
    %[~,pp1]=princom(p1);
    scatter3(p2(:,1),p2(:,2)+90,p2(:,3),5,'r')
    for i=1:post.res.nslice
        plot3(post.res.per(i).points(:,1),post.res.per(i).points(:,2)+90,post.res.per(i).points(:,3),'r','LineWidth',2)
    end
    %scatter3(pp1(:,1),pp1(:,2)+150,pp1(:,3),5,'r')
    axis auto;
    axis equal;
    title(['...' experiment(1).prename(end-20:end)], 'Interpreter', 'none');
end

end

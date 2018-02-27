function [ PTS,T,Pts,Th,score,error,ruPTS] = smooth_periodic( pts,options)
% Finds a smooth periodic line running through a cloud of points
%   pts is a npx3 vector
%   WW is a 1xnp vector containing the weights of each point
%   dt is the step in theta we want to have
%   npp is the number of data points per period we want to use to spline
%   mode and geom determine the fitting procedure. 
%   default mode : spline
%   default geometry : polar
%
% OUTPUT :
% PTS,T : fitted points and polar angle of the smoothed points
% PTS are turned in a new reference frame
% Pts,Th : points and polar angle of experimental points
% score : out of planeness 
% error : standard deviation between smoothed pts and exp pts
%
% Serge Dmitrieff, EMBL, 2015


if nargin<2
    options=pombe_default_options();
end

np=size(pts,1);
Pts=zeros(3,np);
Pts(1:2,:)=pts';
W=ones(1,np);
dt=options.spline.dt;
npp=options.spline.npp;
% --------------------------------------------------------------------
%% Pre-processsing
% --------------------------------------------------------------------
%% Finding mean point and mean plane
% Finding centroid and mean plane of the data.
%[pts,translation]=center_points(Pts,W);
%[ angs,~,pts ] = find_orientation3D( pts,W );
% Converting to cylindrical coordinates
cylpts=cart_to_cyl_z(Pts);
Th=cylpts(1,:);
DPr=cylpts(2,:);
PtZ=cylpts(3,:);
%% Reording the data by increasing angle ; Th is now our coordinate
[Th(1,:),order]=sort(Th);
W(1,:)=W(1,order);
Nr=DPr(1,order); 
NPtZ=PtZ(1,order);
%% Data duplication to enforce periodicity
NW=[W(1,:) W(1,:) W(1,:)];
t=[(-2*pi+Th(1,:)) Th Th(1,:)+2*pi];
% Cartesian coordinates
pts=Pts(:,order); 
% Polar coordinates
r=[Nr(1,:) Nr(1,:) Nr(1,:)];
U=[NPtZ(1,:) NPtZ(1,:) NPtZ(1,:)];
T=-pi:dt:pi;


% --------------------------------------------------------------------
%% Smoothing 
% --------------------------------------------------------------------
% We always convert back to cartesian for simplicity
[ ruPTS ] = spline_fit_periodic( [r;U],NW,t,T,np,npp );
PTS=cyl_to_cart([T;ruPTS]);

% --------------------------------------------------------------------
%% Plotting for debug
% --------------------------------------------------------------------
if 0
	%figure
	%hold all
	%scatter3(turnedpts(1,:),turnedpts(2,:),turnedpts(3,:),'b')
	%plot3(newpts(1,:),newpts(2,:),newpts(3,:),'Linewidth',2)
	figure
	hold all
	%scatter3(NPts(1,:),NPts(2,:),NPts(3,:),'b')
    scatter3(Pts(1,:),Pts(2,:),Pts(3,:),'b')
	%plot3(nPTS(1,:),nPTS(2,:),nPTS(3,:),'Linewidth',2)
    plot3(PTS(1,:),PTS(2,:),PTS(3,:),'Linewidth',2)
end

%cpts=PTS(:,np+1:2*np);
%ppts=newpts(:,np+1:2*np);
end


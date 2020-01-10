function [ PTS ] = spline_fit_periodic( pts,W,t,T,np,npp )
%Fits a smooth periodic line to an order set of points
%   pts is a NDxNnp vector (Nnp: #points, ND : #dimensions)
%   W is a 1xNnp vector containing the weights of each point
%   dt is the step in theta we want to have
%   npp is the number of data points per period we want to use to spline
%   np is the number of points in a period
%
% Serge Dmitrieff, EMBL, 2015

% Preparing variables
WT=0;
NT=length(T);
[ND,Nnp]=size(pts);
% Points we want to interpolate
PTS=zeros(ND,NT);
%number of points per period :
nop=floor(Nnp/3);
dp=floor(nop/npp);

if npp<6
    warning(['Using ' num2str(npp) ' points sampling to spline a period... This is quite few don''t you think ?'])
end

%% Now we're gonna spline and average
% This could be slightly (?) improved by using splinefit
% we're creating splines from about npp points per period
% and then averaging the splines
% based on periodicity, we need to do take around np/npp starting points

npeff=ceil(np/(1*npp));
indexes=1:Nnp;    
for i=1:npeff
    ixes=indexes(i:dp:Nnp);
    %ixes=i:dp:Nnp;
    w=sum(W(ixes));
    for d=1:ND
        PTS(d,:)=PTS(d,:)+spline(t(ixes),pts(d,ixes),T)*w;
    end
    WT=WT+w;
    if 0
        PP(1,:)=spline(t(ixes),pts(1,ixes),T);
        PP(2,:)=spline(t(ixes),pts(2,ixes),T);
        
        figure
        scatter(t(1,:),pts(1,:),10,'b')
        hold all
        scatter(t(1,ixes),pts(1,ixes),80,'r','LineWidth',4)
        plot(T,PP(1,:),'r','LineWidth',3)
        %plot(T-2*pi,PP(1,:),'r','LineWidth',2)
        %plot(T+2*pi,PP(1,:),'r','LineWidth',2)
        figure
        scatter(t(1,:),pts(2,:),10,'b')
        hold all
        scatter(t(1,ixes),pts(2,ixes),80,'r','LineWidth',4)
        plot(T,PP(2,:),'r','LineWidth',3)
        %plot(T-2*pi,PP(2,:),'r','LineWidth',2)
        %plot(T+2*pi,PP(2,:),'r','LineWidth',2)
        
    end
end
%XX=XX/WT;
%YY=YY/WT;
%ZZ=ZZ/WT;
PTS=PTS/WT;


end


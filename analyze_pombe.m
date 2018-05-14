function [ analysis,points] = analyze_pombe( pombe,options)
%analyze_pombe analyzes the shape of pombe cells
%   Analyzes a pombe segmentation
%   Nothing much for now
%
% Serge Dmitrieff, IJM 2018
% www.biophysics.fr

if nargin<2
    options=pombe_default_options();
end

points=pombe.points;
normals=pombe.normals;
faces=pombe.faces;


%% Pre-treatment
% Centering of the data...

if options.centering > 0
    np=size(points,1);
    center=mean(points,1);
    points=points-ones(np,1)*center;
    [~,points,~]=princom(points);
end
%np=size(points,1);
%center=mean(points,1);
%points=points-ones(np,1)*center;

%% We should rotate points & normals
% We use PCA to get the eigen vectors of the points
%[coeffs]=princom(points);
% The eigenvectors define a rotation matrix !
%points=points*coeffs;
%normals=normals*coeffs;
% ...
% Oh wait we don't need the eigen vectors
% It directly gives us the rotated points !
%[coeffs,points,latXYZ]=princom(points);
% We might need to remember the rotation matrix though, when we do comparisons !
%analysis.coeffs=coeffs;

%% We can already compute surface area.
[analysis.surface,analysis.volume]=get_surface_volume(points,faces);

%% Now that it's turned, we can compute one or a couple thin slices
analysis.slices=compute_perimeters(points,options);
%h=options.thickness_central/2.0;
ix=round(analysis.slices.nslice/2);
latYZ=analysis.slices.latXY(ix,:);
%gl=logical((points(:,1)<h).*(points(:,1)>-h));
%slice=points(gl,:);
% We do a PCA of the slice on the YZ plane
%[~,sliceYZ,latYZ]=princom(slice(:,2:3));
% Random definition of circularity
% This is quite sensitive cause we take variance of points
% and then sqrt of ratio -> maximum sensitivity !
analysis.central_circularity=1-sqrt(abs(latYZ(1)-latYZ(2))/sum(latYZ));
%analysis.total_circularity=1-sqrt(abs(latXYZ(3)-latXYZ(3))/sum(latXYZ(2:3)));

analysis.central_r1=sqrt(latYZ(1));
analysis.central_r2=sqrt(latYZ(2));
analysis.total_l1=max(points(:,1))-min(points(:,1));
analysis.pt_mean_dist=mean_dist(points,faces);
% We fit the central slice by a smooth periodic function
%per= smooth_periodic( sliceYZ,options);
%analysis.central_perimeter=seglength(per);
%per=[per,per(:,end)];
analysis.central_perimeter=analysis.slices.pers(ix);
analysis.central_area=analysis.slices.ares(ix);

analysis.length=max(points(:,1))-min(points(:,1));

%% Plotting if we need to
if options.verbose>2
    disp(['Volume : ' num2str(analysis.volume) '   ; surface : ' num2str(analysis.surface)]);
    figure
    hold all
    scatter3(points(:,1),points(:,2),points(:,3),5,'k')
    scatter3(slice(:,1),slice(:,2),slice(:,3),15,'r')
    axis equal
    figure
    scatter(sliceYZ(:,1),sliceYZ(:,2),'r')
    hold all
    plot(per(1,:),per(2,:),'k')
    axis equal
end



end


function [surf]=get_surface_slice(per)
    segs=diff(per,1,2);
    %pp=per(:,1:end-1);
    vecs=cross(segs,per(:,1:end-1),1);
    surf=sum(sqrt(sum(vecs.^2,1)));
    %surf=sum(vecs.^2)/2.0;
end


function [surf,vol]=get_surface_volume(points,face)
    nf=size(face,1);
    surf=0;
    vol=0;
    % Marginal gain
    A=zeros(1,3);
    B=zeros(1,3);
    C=zeros(1,3);
    dir=zeros(1,3);
    pos=zeros(1,3);
    for f=1:nf
        A(:)=points(face(f,1),:);
        B(:)=points(face(f,2),:);
        C(:)=points(face(f,3),:);
        dir(:)=cross((B-A),(C-A));
        surf=surf+norm(dir);
        pos(:)=(A+B+C);
        % This is so freacking cool !!!
        vol=vol+abs(dot(pos,dir));
    end
    surf=surf/2; 
    vol=vol/12.0;
end

function [dist]=mean_dist(points,face)
    nf=size(face,1);
    np=0;
    dist=0;
    % Marginal gain
    A=zeros(1,3);
    B=zeros(1,3);
    C=zeros(1,3);
    for f=1:nf
        A(:)=points(face(f,1),:);
        B(:)=points(face(f,2),:);
        C(:)=points(face(f,3),:);
        dist=dist+sqrt(sum((A-B).^2));
        dist=dist+sqrt(sum((C-B).^2));
        dist=dist+sqrt(sum((A-C).^2));
        np=np+3;
    end
    dist=dist/np;
end


function [res]=compute_perimeters(points,options)
    h=options.thickness_central/2.0;
    mp=min(points(:,1));
    Mp=max(points(:,1));
    N=floor((Mp-mp-2*h)/(2*h));
    res.pers=zeros(N,1);
    res.ares=zeros(N,1);
    res.pos=zeros(N,1);
    res.latXY=zeros(N,2);
    res.per(N).points=[0 0 0];
    xm=mp+2*h;
    %figure
    %hold all
    for i=1:N
        res.pos(i)=xm+((i-1)*h*2);
        
        gl=logical((points(:,1)<(xm+(i*h*2))).*(points(:,1)>(xm-2*h+(i*h*2))));
        slice=points(gl,:);
        [~,~,res.latYZ(i,:)]=princom(slice(:,2:3));
        %scatter3(slice(:,1),slice(:,2),slice(:,3))
        per= smooth_periodic( slice(:,2:3),options);
        %
        per=[per,per(:,1)];
        
        res.ares(i)=get_surface_slice(per);
        res.pers(i)=seglength(per);
        res.per(i).points=[res.pos(i)*ones(size(per,2),1) per(1:2,:)'];
    end
    res.nslice=N;
end
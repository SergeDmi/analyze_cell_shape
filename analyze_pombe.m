function [ analysis] = analyze_pombe( pombe,options)
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
np=size(points,1);
center=mean(points,1);
points=points-ones(np,1)*center;

%% We should rotate points & normals
% We use PCA to get the eigen vectors of the points
%[coeffs]=princom(points);
% The eigenvectors define a rotation matrix !
%points=points*coeffs;
%normals=normals*coeffs;
% ...
% Oh wait we don't need the eigen vectors
% It directly gives us the rotated points !
[~,points]=princom(points);
% We might need to remember the rotation matrix though, when we do comparisons !


%% We can already compute surface area.
[analysis.surface,analysis.volume]=get_surface_volume(points,faces);

%% Now that it's turned, we can compute one or a couple thin slices
h=options.thickness_central/2.0;
gl=logical((points(:,1)<h).*(points(:,1)>-h));
slice=points(gl,:);
% We do a PCA of the slice on the YZ plane
[~,sliceYZ,latYZ]=princom(slice(:,2:3));
% Random definition of circularity
% This is quite sensitive cause we take variance of points
% and then sqrt of ratio -> maximum sensitivity !
analysis.circularity=1-sqrt(abs(latYZ(1)-latYZ(2))/sum(latYZ));
analysis.central_r1=sqrt(latYZ(1));
analysis.central_r2=sqrt(latYZ(2));
analysis.pt_mean_dist=mean_dist(points,faces);
% We fit the central slice by a smooth periodic function
per= smooth_periodic( sliceYZ,options);
analysis.perimeter=seglength(per);
per=[per,per(:,end)];
analysis.slice_area=get_surface_slice(per);

analysis.length=max(points(:,1))-min(points(:,1));

%% Plotting if we need to
if options.verbose>0
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


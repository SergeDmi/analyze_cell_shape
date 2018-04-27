function [ analysis] = integrate_pombe( pombe,options)
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
R0=options.R0;

%% Pre-treatment
% Centering of the data...
np=size(points,1);
nf=size(faces,1);
center=mean(points,1);
points=points-ones(np,1)*center;
[coefs,points,~]=princom(points);
normals=normals*coefs;

% Curvature along s
Cs=zeros(np,1);
% Curvature along phi
Cp=zeros(np,1);

%recomputing normals
norms=get_normals(points,faces,np,nf);
% forces
Fpts=zeros(np,3);
% Interaction matrix
matR=create_interaction_matrix(points,norms,faces);


%% We should rotate points & normals
% We use PCA to get the eigen vectors of the points
%[coeffs]=princom(points);
% The eigenvectors define a rotation matrix !
%points=points*coeffs;
%normals=normals*coeffs;
% ...
% Oh wait we don't need the eigen vectors
% It directly gives us the rotated points !

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
analysis.central_circularity=1-sqrt(abs(latYZ(1)-latYZ(2))/sum(latYZ));
analysis.total_circularity=1-sqrt(abs(latXYZ(3)-latXYZ(3))/sum(latXYZ(2:3)));

analysis.central_r1=sqrt(latYZ(1));
analysis.central_r2=sqrt(latYZ(2));
analysis.total_r1=sqrt(latXYZ(2));
analysis.total_r2=sqrt(latXYZ(3));
analysis.total_l1=sqrt(latXYZ(1));
analysis.pt_mean_dist=mean_dist(points,faces);
% We fit the central slice by a smooth periodic function
per= smooth_periodic( sliceYZ,options);
analysis.central_perimeter=seglength(per);
per=[per,per(:,end)];
analysis.central_area=get_surface_slice(per);

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


function [M]=normalize_rows_3D(M)
M(:,:)=M./(sqrt(sum(M.^2,2))*[1 1 1]);
end

function [matR]=create_interaction_matrix(points,faces,normsP,normsF)
% initiation
np=size(points,1);
nf=size(faces,1);
matR=zeros(np,np);
dir=ones(np,1)*[0 0 1];
% we create the basis
n_psi=normalize_rows_3D(cross(normals,dir,2));
n_sss=cross(n_psi,n_sss,2);

CP=zeros(nf,3);
CS=zeros(nf,3);
% Now we need to compute curvatures !
% Marginal gain
A=zeros(1,3);
B=zeros(1,3);
C=zeros(1,3);
dir=zeros(1,3);
pos=zeros(1,3);

gradN=zeros(3,3);
NP=zeros(3,1);
NS=zeros(3,1);
for f=1:nf
	iA=face(f,1);
	iB=face(f,2);
	iC=face(f,3);
	A(:)=points(iA,:);
	B(:)=points(iB,:);
	C(:)=points(iC,:);
	NP=
	%% 
	gradN(:,:)=surface_grad(A,B,normsP(iA,:),normsP(iB,:));
	C_p=gradN(n_psi(iA,:)+n_psi(iB,:))'/2.0;
	C_p=(n_psi(iA,:)+n_psi(iB,:))'/2.0;
	
	dir(:)=cross((B-A),(C-A));
	surf=surf+norm(dir);
	pos(:)=(A+B+C);
	% This is so freacking cool !!!
	vol=vol+abs(dot(pos,dir));
end
	

end

function dNdAB=surface_grad(A,B,nA,nB)
%invdir=1./(B-A);
%dN=nB-nA;
dNdAB=(1./(B-A)')*(nB-nA);
end

function [normsP,normF]=get_normals(points,faces,np,nf)
    %nf=size(face,1);
    %surf=0;
	%vol=0;
    % Marginal gain
	normsP=zeros(np,3);
	normsF=zeros(np,3);
    %A=zeros(1,3);
    %B=zeros(1,3);
    %C=zeros(1,3);
    dir=zeros(1,3);
	col3=ones(3,1);
    for f=1:nf
        %A(:)=points(faces(f,1),:);
        %B(:)=points(faces(f,2),:);
        %C(:)=points(faces(f,3),:);
        %dir(:)=cross((B-A),(C-A));
		dir(:)=cross((points(faces(f,2),:)-points(faces(f,1),:)),(points(faces(f,3),:)-points(faces(f,1),:)));
		normsP(faces(f),:)=normsP(faces(f),:)+col3*dir;
		normsF(f,:)=dir;
    end
    normsP(:)=normalize_rows_3D(normsF);
	normsF(:)=normalize_rows_3D(normsP);
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


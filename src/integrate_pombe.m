function [ points,matR,Cp,Cs,links,norms,ono] = integrate_pombe( pombe,options)
%analyze_pombe analyzes the shape of pombe cells
%   Analyzes a pombe segmentation
%   Nothing much for now
%
% Serge Dmitrieff, IJM 2018
% www.biophysics.fr


%% Checking inputs
if nargin<2
    options=pombe_default_options();
end

points=pombe.points;
normals=pombe.normals;
faces=pombe.faces;
R0=3;

%% Pre-treatment
% Centering of the data...
np=size(points,1);
nf=size(faces,1);
center=mean(points,1);
points=points-ones(np,1)*center;
%[coefs,points,~]=get_pca_cov(points,1);
[coefs,points,~]=princom(points);
normals=normals*coefs;

%recomputing normals
[norms,surfs]=get_normals_surfs(points,faces,np,nf);
% forces
Fpts=zeros(np,3);
% Interaction matrix
[matR]=create_interaction_matrix(points,norms,faces);







end

%% Normalize a column of row vectors
function [M,w]=normalize_rows_3D(M)
w=sqrt(sum(M.^2,2));
M(:,:)=M./(w*[1 1 1]);
end

%% Create the interaction matrix
function [matR]=create_interaction_matrix(points,normsP,faces,params)
% initiation
np=size(points,1);
nf=size(faces,1);
matR=zeros(np,np);
dir=ones(np,1)*[1 0 0];
% we create the basis
n_psi=normalize_rows_3D(cross(normsP,dir,2));
%n_sss=cross(n_psi,normsP,2);
% Marginal gain
A=zeros(1,3);
B=zeros(1,3);
NP=zeros(3,1);
NS=zeros(3,1);





for f=1:nf

	for i=1:3
		if i==1
			iA=faces(f,1);
			iB=faces(f,2);
		elseif i==2
			iA=faces(f,1);
			iB=faces(f,3);
		elseif i==3
			iA=faces(f,2);
			iB=faces(f,3);
		end

        if matR(iA,iB)==0 && matR(iB,iA)==0
           matA(iA,iB)=create_link(iA,iB,points,n_psi,params);
        end
    end
end


end

function l0=create_link(iA,iB,points,n_psi,params)
% Dummy
BA(:)=points(iB,:)-points(iA,:);
%B(:)=points(iB,:);
NP(:)=(n_psi(iA,:)+n_psi(iB,:))'/2.0;

l0=norm(BA)
%NS(:)=(n_sss(iA,:)+n_sss(iB,:))'/2.0;



end

function [normsP]=get_normals(points,faces,np,nf)
	normsP=zeros(np,3);
	normsF=zeros(np,3);
    dir=zeros(1,3);
	col3=ones(3,1);
    for f=1:nf
		dir(:)=cross((points(faces(f,2),:)-points(faces(f,1),:)),(points(faces(f,3),:)-points(faces(f,1),:)));
        normsP(faces(f,:),:)=normsP(faces(f,:),:)+col3*dir;
    end
    normsP(:,:)=normalize_rows_3D(normsP);

end


function [normsP,surfP]=get_normals_surf(points,faces,np,nf)
    % Marginal gain

	normsP=zeros(np,3);
	cntfP=zeros(np,1);
    dir=zeros(1,3);
	col3=ones(3,1);


    for f=1:nf
		dir(:)=cross((points(faces(f,2),:)-points(faces(f,1),:)),(points(faces(f,3),:)-points(faces(f,1),:)));

		normsP(faces(f,:),:)=normsP(faces(f,:),:)+col3*dir;
		cntfP(faces(f,:))=cntfP(faces(f,:))+1;

    end
    [normsP(:,:),surfP]=normalize_rows_3D(normsP);
	surfP(:)=surfP./cntfP;

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

function dummy_function()
KR=(zeros(3*np,3*np));
MR=(zeros(3*np,3*np));
MPP=(zeros(3*np,3*np));
DP=(zeros(3*np,3*np));
DP2=(zeros(3*np,3*np));
FDP=(zeros(3*np,3*np));
PR=zeros(3*np,1);
NDP=(zeros(3*np,3*np));

DEF=ones(3*np,3*np);
NR=zeros(3*np,1);
COL=ones(3*np,1);
ROW=ones(1,3*np);
I=eye(3);
t=0;
ixes=reshape(1:(np*np*9),3*np,3*np);
for i=0:2
	MR((1:np)*3+i-2,(1:np)*3+i-2)=matR(:,:);
	KK=zeros(3*np,3*np);
	KK((1:np)*3+i-2,(1:np)*3+i-2)=matR(:,:)>0;
	if i==0
		allX=ixes(logical(KK));
	elseif i==1
		%KR((1:np)*3+i-2,(1:np)*3+i-2)=matR(:,:);
		KR(:,:)=KK(:,:);
		allY=ixes(logical(KK));
	elseif i==2
		allZ=ixes(logical(KK));
	end
end

PR(:)=reshape(points',[3*np,1]);
OR=MR>0;
KR(:)=KR>0;

%

DEF=DEF-OR;

dt=0.01;
P=1.0;
L0=3;
Tend=0;
Tnorm=10*dt;
nextN=0;
%% To do for all times
while t<Tend
	t=t+dt
	nextN=nextN-dt;
	if nextN<0
		points(:,:)=reshape(PR,[3,np])';
		disp('finding normals')
		tic
		NR(:)=reshape(get_normals(points,faces,np,nf)',[3*np,1]);
		toc
		nextN=Tnorm;
	end
	disp('main loop')
	tic
	MPP(:,:)=OR.*(PR*ROW);
	DP(:,:)=transpose(MPP)-MPP;
	DP2(:,:)=DP.^2;
	NDP(:,:)=DEF;
	NDP(allY)=DP2(allY)+DP2(allX)+DP2(allZ);
	NDP(allX)=NDP(allY);
	NDP(allZ)=NDP(allY);
	FDP(:,:)=(DP-L0*DP./NDP);
	PR(:)=PR(:)+dt*(sum(FDP,2)+NR*P);
	toc

end
end

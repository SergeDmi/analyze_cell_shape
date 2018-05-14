function [ points,matR,Cp,Cs,links,norms,ono] = integrate_pombe( pombe,options)
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
R0=3;

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
[matR,Cp,Cs,links,norms,ono]=create_interaction_matrix(points,norms,faces);

KR=zeros(3*np,3*np);
MR=zeros(3*np,3*np);
MPP=zeros(3*np,3*np);
DP=zeros(3*np,3*np);
DP2=zeros(3*np,3*np);
FDP=zeros(3*np,3*np);
PR=zeros(3*np,1);
NDP=zeros(3*np,3*np);
DEF=ones(3*np,3*np);
NR=zeros(3*np,1);
COL=ones(3*np,1);
ROW=ones(1,3*np);
I=eye(3);
t=0;
for i=0:2
	MR((1:np)*3+i-2,(1:np)*3+i-2)=matR(:,:);
	if i==1
		KR((1:np)*3+i-2,(1:np)*3+i-2)=matR(:,:);
	end
end

PR(:)=reshape(points',[3*np,1]);
OR=MR>0;
KR(:)=KR>0;
DEF=DEF-OR;

dt=0.01;
P=1.0;
L0=3;
Tend=1;
%% To do for all times
while t<Tend
	t=t+dt
	NR(:)=reshape(get_normals(points,faces,np,nf)',[3*np,1]);
	MPP(:,:)=OR.*(PR*ROW);
	DP(:,:)=transpose(MPP)-MPP;
	DP2(:,:)=DP.^2;
	NDP(:,:)=DEF+conv2(sqrt(conv2(DP2,I,'same').*KR),I,'same');
	FDP(:,:)=(DP-L0*DP./NDP);
	PR(:)=PR(:)+dt*(sum(FDP,2)+NR*P);
	points(:,:)=reshape(PR,[3,np])';
end

end


function [M]=normalize_rows_3D(M)
M(:,:)=M./(sqrt(sum(M.^2,2))*[1 1 1]);
end

function [matR,CP,CS,links,norms,ono]=create_interaction_matrix(points,normsP,faces)
% initiation
np=size(points,1);
nf=size(faces,1);
CP=zeros(np,1);
CS=zeros(np,1);
wCP=zeros(np,1);
wCS=zeros(np,1);
matR=zeros(np,np);
dir=ones(np,1)*[1 0 0];
links=zeros(3*nf,7);
count_l=0;
% we create the basis
n_psi=normalize_rows_3D(cross(normsP,dir,2));
n_sss=cross(n_psi,normsP,2);

% Now we need to compute curvatures !
% Marginal gain
A=zeros(1,3);
B=zeros(1,3);

gradN=zeros(3,3);
NP=zeros(3,1);
NS=zeros(3,1);

norms=zeros(3*nf,1);
ono=zeros(3*nf,1);
for f=1:nf
	%iA=face(f,1);
	%iB=face(f,2);
	%iC=face(f,3);
	%A(:)=points(iA,:);
	%B(:)=points(iB,:);
	%C(:)=points(iC,:);
	%NP=n_psi(iA,:)'+n_psi(iB,:)'+n_psi(iC,:)';
	%NS=n_sss(iA,:)'+n_sss(iB,:)'+n_sss(iC,:)';
	
	% Link A to B
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
		A(:)=points(iA,:);
		B(:)=points(iB,:);
		gradN(:,:)=surface_grad(A,B,normsP(iA,:),normsP(iB,:));
		NP(:)=(n_psi(iA,:)+n_psi(iB,:))'/2.0;
		NS(:)=(n_sss(iA,:)+n_sss(iB,:))'/2.0;
        % Weighted average of curvatures
		CP(iA)=CP(iA)+dot(gradN*NP,NP)*(dot(NP,B-A)^2);
		CS(iA)=CS(iA)+dot(gradN*NS,NS)*(dot(NS,B-A)^2);
        CP(iB)=CP(iB)+dot(gradN*NP,NP)*(dot(NP,B-A)^2);
		CS(iB)=CS(iB)+dot(gradN*NS,NS)*(dot(NS,B-A)^2);
        wCP(iA)=wCP(iA)+(dot(NP,B-A)^2);
        wCP(iB)=wCP(iB)+(dot(NP,B-A)^2);
        wCS(iA)=wCS(iA)+(dot(NS,B-A)^2);
        wCS(iB)=wCS(iB)+(dot(NS,B-A)^2);
		%Ss=1.0/(2.0*Cp);
		%Sp=(1.0-Cs*Ss)/Cp;
		%matR(iA,iB)=(Ss*dot(B-A,NS)+Sp*dot(B-A,NP))/norm(B-A);
		%CP(iA)=Cp;
		%CS(iA)=Cs;
		
	end
end
CP=CP./wCP;
CS=CS./wCS;

for f=1:nf
	%iA=face(f,1);
	%iB=face(f,2);
	%iC=face(f,3);
	%A(:)=points(iA,:);
	%B(:)=points(iB,:);
	%C(:)=points(iC,:);
	%NP=n_psi(iA,:)'+n_psi(iB,:)'+n_psi(iC,:)';
	%NS=n_sss(iA,:)'+n_sss(iB,:)'+n_sss(iC,:)';
	
	% Link A to B
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
		A(:)=points(iA,:);
		B(:)=points(iB,:);
        NP(:)=(n_psi(iA,:)+n_psi(iB,:))'/2.0;
		NS(:)=(n_sss(iA,:)+n_sss(iB,:))'/2.0;
        Ss=1.0/(CP(iA)+CP(iB));
		Sp=2.0*(1.0-0.5*(CS(iA)+CS(iB))*Ss)/(CP(iA)+CP(iB));
        if matR(iA,iB)==0 && matR(iB,iA)==0
            %matR(iA,iB)=(Ss*dot(B-A,NS)+Sp*dot(B-A,NP))/norm(B-A);
            %matR(iB,iA)=matR(iA,iB);
            count_l=count_l+1;
            %size([A(:)',B(:)',matR(iA,iB)])
            %links(count_l,:)=[A(:)',B(:)',matR(iA,iB)];
			%norms(count_l)=norm(B-A);
			%ono(count_l)=(dot(B-A,NP)/norms(count_l)).^2;
			matR(iA,iB)=1;
            matR(iB,iA)=1;
        end
    end
end

links=links(1:count_l,:);
ono=ono(1:count_l,:);
norms=norms(1:count_l,:);

if 0
    figure
    hold all
    scatter3(points(1:300,1),points(1:300,2),points(1:300,3),'k')
    for n=1:300
        plot3([points(n,1) points(n,1)-5*normsP(n,1)],[points(n,2) points(n,2)-5*normsP(n,2)],[points(n,3) points(n,3)-5*normsP(n,3)],'k');
        plot3([points(n,1) points(n,1)-5*n_psi(n,1)],[points(n,2) points(n,2)-5*n_psi(n,2)],[points(n,3) points(n,3)-5*n_psi(n,3)],'r');
        plot3([points(n,1) points(n,1)-5*n_sss(n,1)],[points(n,2) points(n,2)-5*n_sss(n,2)],[points(n,3) points(n,3)-5*n_sss(n,3)],'g');
    end	
end

end

function dNdAB=surface_grad(A,B,nA,nB)
%invdir=1./(B-A);
%dN=nB-nA;
dNdAB=(1./(B-A)')*(nB-nA);
end

function [normsP,normsF]=get_normals(points,faces,np,nf)
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
        
		normsP(faces(f,:),:)=normsP(faces(f,:),:)+col3*dir;
		normsF(f,:)=dir;
    end
    normsP(:,:)=normalize_rows_3D(normsP);
	normsF(:,:)=normalize_rows_3D(normsF);
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


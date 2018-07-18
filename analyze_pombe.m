function [ analysis,points] = analyze_pombe( pombe,options)
%analyze_pombe analyzes the shape of pombe cells
%   Analyzes a pombe segmentation
%   Nothing much for now
%
% Serge Dmitrieff, IJM 2018
% www.biophysics.fr

if nargin<2
    options=pombe_default_options();
    if isfield(options,"pixel_size")
      pixel_size=options.pixel_size;
    else
      disp('Warning : could not find pixel_size in options');
      pixel_size=1;
    end
end

points=options.pixel_size*pombe.points;
normals=pombe.normals;
faces=pombe.faces;


%% Pre-treatment
% Centering of the data...
np=size(points,1);
nf=size(faces,1);
if options.centering > 0
    center=mean(points,1);
    points=points-ones(np,1)*center;
    [~,points,~]=princom(points);
	   disp('Auto centering and aligning cell')
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


[normsP,normsF,surfP,surfF,matR]=get_normals(points,faces,np,nf);
[CP,CS]=compute_curvatures(points,faces,normsP,np,nf);
analysis.surface_points=sum(surfP);
analysis.surface_faces=sum(surfF);

momsS=moments(CS,3);
momsP=moments(CP,3);
[ofs,ks,ths]=gamma_params(momsS);
[ofp,kp,thp]=gamma_params(momsP);
analysis.distrCS=[ofs,ks,ths];
analysis.distrCP=[ofp,kp,thp];
analysis.CS=CS;
analysis.CP=CP;
analysis.mean_curv=CP+CS;
analysis.curv_energy=sum(exclude_extremes((analysis.mean_curv.^2).*surfP,3));


%% Now that it's turned, we can compute one or a couple thin slices
analysis.slices=compute_perimeters(points,options);
%h=options.thickness_central/2.0;
ix=round(analysis.slices.nslice/2);
latYZ=analysis.slices.latYZ(ix,:);
%gl=logical((points(:,1)<h).*(points(:,1)>-h));
%slice=points(gl,:);
% We do a PCA of the slice on the YZ plane
%[~,sliceYZ,latYZ]=princom(slice(:,2:3));
% Random definition of circularity
% This is quite sensitive cause we take variance of points
% and then sqrt of ratio -> maximum sensitivity !
analysis.central_circularity=1-sqrt(abs(latYZ(1)-latYZ(2))/sum(latYZ));
%analysis.total_circularity=1-sqrt(abs(latXYZ(3)-latXYZ(3))/sum(latXYZ(2:3)));
analysis.circularity_fiji=1-(sqrt(1-(4*pi*(get_surface_slice(per)/(seglength(per)).^2))));
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

function [M,w]=normalize_rows_3D(M)
w=sqrt(sum(M.^2,2));
M(:,:)=M./(w*[1 1 1]);
end

function [surf]=get_surface_slice(per)
    segs=diff(per,1,2);
    %pp=per(:,1:end-1);
    vecs=cross(segs,per(:,1:end-1),1);
    surf=sum(sqrt(sum(vecs.^2,1)));
    %surf=sum(vecs.^2)/2.0;
end

function curvs=get_mean_curv(points,normsP,matR,np)
curvs=zeros(np,1);
coefs=zeros(np,1);
for i=2:np
	inters=logical(matR(i,1:(i-1)));
	ni=sum(inters);
	dist2s=sum((ones(ni,1)*points(i,:)-points(inters,:)).^2,2);
	prods=(1-sum((ones(ni,1)*normsP(i,:)).*points(inters,:),2))./dist2s;
	curvs(inters)=curvs(inters)+sqrt(prods);
	coefs(inters)=coefs(inters)+1;
	curvs(i)=curvs(i)+sum(prods);
	coefs(i)=coefs(i)+1;
end
curvs=curvs./coefs;
end

function dNdAB=surface_grad(A,B,nA,nB)
%invdir=1./(B-A);
%dN=nB-nA;
dNdAB=(1./(B-A)')*(nB-nA);
end

function [normsP,normsF,surfP,surfF,matR]=get_normals(points,faces,np,nf)
    % Marginal gain
	matR=zeros(np,np);
	normsP=zeros(np,3);
	normsF=zeros(np,3);
	%surfF=zeros(nf,3);
	%surfP=zeros(np,3);
	cntfP=zeros(np,1);
	cntfF=zeros(nf,1);
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
		matR(faces(f,1),faces(f,2))=1;
		matR(faces(f,2),faces(f,3))=1;
		matR(faces(f,3),faces(f,1))=1;
        %surfP(faces(f,:),:)=surfP(faces(f,:),:)+
		normsP(faces(f,:),:)=normsP(faces(f,:),:)+col3*dir;
		cntfP(faces(f,:))=cntfP(faces(f,:))+1;
		normsF(f,:)=dir;
		cntfF(f)=cntfF(f)+1;


		% Link A to B


    end
    [normsP(:,:),surfP]=normalize_rows_3D(normsP);
	[normsF(:,:),surfF]=normalize_rows_3D(normsF);
	surfP(:)=surfP./cntfP;
	surfF(:)=surfF./cntfF;


end

function [CP,CS]=compute_curvatures(points,faces,normsP,np,nf)
poins(np).lCs=[];

gradN=zeros(3,3);
NP=zeros(3,1);
NS=zeros(3,1);

dirx=ones(np,1)*[1 0 0];
% we create the basis
n_psi=normalize_rows_3D(cross(normsP,dirx,2));
n_sss=cross(n_psi,normsP,2);

for i=1:np
	poins(i).lCs=[];
	poins(i).lCp=[];
	poins(i).wCs=[];
	poins(i).wCp=[];
end
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
		A(:)=points(iA,:);
		B(:)=points(iB,:);
		gradN(:,:)=surface_grad(A,B,normsP(iA,:),normsP(iB,:));
		NP(:)=(n_psi(iA,:)+n_psi(iB,:))'/2.0;
		NS(:)=(n_sss(iA,:)+n_sss(iB,:))'/2.0;
		% Weighted average of curvatures

		wp=(dot(NP,B-A)^2);
		ws=(dot(NS,B-A)^2);
		poins(iA).lCp=[poins(iA).lCp dot(gradN*NP,NP)];
		poins(iA).lCs=[poins(iA).lCs dot(gradN*NS,NS)];
		poins(iB).lCp=[poins(iB).lCp dot(gradN*NP,NP)];
		poins(iB).lCs=[poins(iB).lCs dot(gradN*NS,NS)];
		poins(iA).wCp=[poins(iA).wCp wp];
		poins(iA).wCs=[poins(iA).wCs ws];
		poins(iB).wCp=[poins(iB).wCp wp];
		poins(iB).wCs=[poins(iB).wCs ws];


	end
end

[CP,CS]=get_curvs(poins,np);
end

function [Cp,Cs]=get_curvs(poins,np)
Cp=zeros(np,1);
Cs=zeros(np,1);
for i=1:np
	[lCp,icp]=unique(poins(i).lCp);
	wCp=poins(i).wCp(icp);
	lCs=poins(i).lCs(icp);
	wCs=poins(i).wCs(icp);

	[lCp,ix]=exclude_extremes(lCp,2);
	wCp=wCp(ix);
	[lCs,ix]=exclude_extremes(lCs,2);
	wCs=wCs(ix);
	Cp(i)=sum(lCp.*wCp)/sum(wCp);
	Cs(i)=sum(lCs.*wCs)/sum(wCs);
end
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
	res.dist=zeros(N,1);
    res.latYZ=zeros(N,2);
    res.per(N).points=[0 0 0];
    xm=mp+2*h;
    %figure
    %hold all
    for i=1:N
		pos=xm+(i*h*2)-h;
        res.pos(i)=pos;
        res.dist(i)=min(abs(mp-pos),abs(Mp-pos));

        gl=logical((points(:,1)<(pos+h)).*(points(:,1)>(pos-h)));
        slice=points(gl,:);
        [~,~,res.latYZ(i,:)]=princom(slice(:,2:3));
        %scatter3(slice(:,1),slice(:,2),slice(:,3))
        per= smooth_periodic( slice(:,2:3),options);

        %
        per=[per,per(:,1)];
        res.circ(i)=1-sqrt(abs(res.latYZ(i,1)-res.latYZ(i,2))/sum(res.latYZ(i,:)));
        res.ares(i)=get_surface_slice(per);
        res.pers(i)=seglength(per);
        res.per(i).points=[res.pos(i)*ones(size(per,2),1) per(1:2,:)'];
    end
    res.nslice=N;
end

function [ analysis,points] = analyze_shape( shape,options)
% analyze_shape analyzes the shape of a cell
%   Computes a set of properties such as volume, surface, length, etc.
%   shape is a container for points (coordinates), faces and normals
%
% Serge Dmitrieff, IJM 2018
% www.biophysics.fr



%% We  first test the input
if nargin<2
    options=analysis_default_options();
end

if isfield(options,'pixel_size')
  pixel_size=options.pixel_size;
else
  disp_if_verbose('Warning : could not find pixel_size in options',options);
  pixel_size=1;
end

%% We resize the points if needed
if options.do_shape_resize
  points=pixel_size*shape.points;
else
  points=shape.points;
end
normals=shape.normals;
faces=shape.faces;


%% Pre-treatment
% Centering of the points...
np=size(points,1);
nf=size(faces,1);
if options.centering > 0
  center=mean(points,1);
  points=points-ones(np,1)*center;
	disp_if_verbose('Auto centering cell',options)
end
% Aligning the point
if options.aligning > 0
	 [~,points,~]=princom(points);
	  disp_if_verbose('Auto-aligning cell',options)
end
if options.do_backbone_analysis
  % Now we're going to draw a central axis through the cell
  [analysis.backbone_points,analysis.backbone_fits,analysis.backbone_arclen]=get_rolling_backbone(points,options.backbone);
  analysis.backbone_curv  =get_backbone_curvs(analysis.backbone_points);
end


%% We can already compute surface area.
[analysis.surface,analysis.volume]=get_surface_volume(points,faces);
[normsP,normsF,surfP,surfF,matR]=get_normals(points,faces,np,nf);

%% We can compute the curvature of the cell wall on the main axes
% CS corresponds to curvature on the cell main axis
% CP corresponds to curvature on the cell short axis
if (options.do_curvature_analysis)
  [CP,CS]=compute_curvatures(points,faces,normsP,np,nf);
  analysis.surface_points=sum(surfP);
  analysis.surface_faces=sum(surfF);

  %momsS=moments(CS,3);
  %momsP=moments(CP,3);
  %[ofs,ks,ths]=gamma_params(momsS);
  %[ofp,kp,thp]=gamma_params(momsP);
  %analysis.distrCS=[ofs,ks,ths];
  %analysis.distrCP=[ofp,kp,thp];
  analysis.CS=CS;
  analysis.CP=CP;
  analysis.mean_curv=CP+CS;
  analysis.curv_energy=sum(exclude_extremes((analysis.mean_curv.^2).*surfP,3));
end

%% We can also chop the cell into think slices and compute their shapes
if (options.do_slice_analysis)
  analysis.slices=make_analyze_slices(points,options);
  %h=options.thickness_central/2.0;
  ix=round(analysis.slices.nslice/2);
  latYZ=analysis.slices.latYZ(ix,:);
  sliceYZ=analysis.slices.per(ix).points(:,2:3);
  %analysis.mean_circ=mean(analysis.slices.circ);
  % Weighted mean
  analysis.mean_circ=sum(analysis.slices.circ.*analysis.slices.ares')/sum(analysis.slices.ares);
  analysis.mean_pers=mean(analysis.slices.pers);
  analysis.mean_area=mean(analysis.slices.ares);
else
  % We need at least the central slice
  h=options.thickness/2.0;
  gl=logical((points(:,1)<h).*(points(:,1)>-h));
  slice=points(gl,:);
  [~,sliceYZ,latYZ]=princom(slice(:,2:3));
end

%% We can approximate the centeral slice using a custom spline fitting method
per=smooth_periodic(sliceYZ,options.central.spline);
per=[per(:,end) per];

%% From this we compute shapes
[perimeter,area,excess_length]=compute_slice_shape(per);
analysis.central_perimeter=perimeter;
analysis.central_area=area;
analysis.perimeter_angular_deviation=compute_angular_deviations(per(1:2,2:end)');
analysis.central_spline=per';
analysis.central_points=sliceYZ;
% This is a bad way to compute circularity....
analysis.central_circularity=1-sqrt(abs(latYZ(1)-latYZ(2))/sum(latYZ));
% This is a better way
analysis.circularity_fiji=1-(sqrt(1-(4*pi*(area/(perimeter^2)))));
% Equivalently, we can compute the excess perimeter compared to a circle of similar area
analysis.excess_length=excess_length;
analysis.central_r1=sqrt(latYZ(1));
analysis.central_r2=sqrt(latYZ(2));
analysis.total_l1=max(points(:,1))-min(points(:,1));
analysis.pt_mean_dist=mean_dist(points,faces);
analysis.length=max(points(:,1))-min(points(:,1));

%@TODO : unified way to compute circularity and straightness...
% Tough ! This is implemented in the compute_angular deviation part, works more or less


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

%% UTILITY FUNCTIONS

% Sort a bunch of points by their angle
function [PTS]=sort_by_angle(pts,options)
  % Assumes PTS as a column vector of points (row of coordinates)
np=size(pts,1);
Pts=zeros(3,np);
Pts(1:2,:)=pts';
% converts to cylindrical coordinates to find the order
cylpts=cart_to_cyl_z(Pts);
[~,order]=sort(cylpts(1,:));
% gotcha
PTS=Pts(1:2,order)';
end

% This function computes the std deviation of the anges between consecutive segments in a set of points
function [devi]=compute_angular_deviations(pts)
% first find the orientation of segments
dp=normalize_rows_ND(diff(pts));
ns=size(dp,1);
% we use the dot product to compute angle
projs=dp(1:ns-1,:).*dp(2:ns,:);
angs=acos(sum(projs,2));
devi=std(angs);
end

% Normalize the rows in a column vector of row vectors
function [M,w]=normalize_rows_3D(M)
w=sqrt(sum(M.^2,2));
M(:,:)=M./(w*[1 1 1]);
end

% We compute the shape of a slice of cell wall
% The slice should be already splined !
function [perimeter,area,excess_length,circ]=compute_slice_shape(per)
perimeter=seglength(per);
area=get_surface_slice(per);
% this is the relative excess perimeter of the slice compared to the circle of same surface area
% with this definition, e_l gives dr/r for an ellipse of radii r,r+dr with dr -> 0
excess_length=4*sqrt((perimeter^2)/(4*pi*area)-1);
circ=1-(sqrt(1-(4*pi*(area/(perimeter^2)))));
% Fiji circularity
end

% Computes the surface area of a smoothed slice
function [surf]=get_surface_slice(per)
% the surface area is the sum of the areas of the triangles between two consecutive points and the center
% i.e. the sum of the cross products !
segs=diff(per,1,2);
vecs=cross(segs,per(:,1:end-1),1);
surf=sum(sqrt(sum(vecs.^2,1)))/2.0;
end

% Compute the surface Jacobian
function dNdAB=surface_grad(A,B,nA,nB)
dNdAB=(1./(B-A)')*(nB-nA);
end

% computes the main curvatures of a discretized curve
function [curvs]=get_backbone_curvs(ske)
% we compute the arclength along the curve
arclen=get_arclength(ske)';
% we align the curve
[~,ske]=princom(ske);
% we compute the curvature of y,z as a function of the arclength
py=polyfit(arclen,ske(:,2),2);
pz=polyfit(arclen,ske(:,3),2);
curvs=[py(1),pz(1)];
end


%  computes arclength along a discretized curve
function [lens]=get_arclength(M)
  lens=get_cumulative_seglength(M);
  lens=[0;lens]';
end

% computes the cumulative segment length of a discretized curve
function [clen]=get_cumulative_seglength(M)
  lens=get_seglenth(M);
  clen=lens;
  for i=2:numel(clen)
    clen(i)=clen(i)+clen(i-1);
  end
end

% computes the length of segments of a discretized curve
function [lens]=get_seglenth(M)
  sm=size(M);
  segs=diff(M,1,1);
  lens=sqrt(sum(segs.^2,2));
end

% Compute the normals to a suface
function [normsP,normsF,surfP,surfF,matR]=get_normals(points,faces,np,nf)
  % matR is an interaction matrix between points
	matR=zeros(np,np);
  % normsP are the normals to points
	normsP=zeros(np,3);
  % normsF are the normals to faces
	normsF=zeros(np,3);

  % cntfP,F are the counts
	cntfP=zeros(np,1);
	cntfF=zeros(nf,1);
  dir=zeros(1,3);
	col3=ones(3,1);


  for f=1:nf
    % the director vector to a face is the cross products of the vectors
		dir(:)=cross((points(faces(f,2),:)-points(faces(f,1),:)),(points(faces(f,3),:)-points(faces(f,1),:)));
		matR(faces(f,1),faces(f,2))=1;
		matR(faces(f,2),faces(f,3))=1;
		matR(faces(f,3),faces(f,1))=1;
		normsP(faces(f,:),:)=normsP(faces(f,:),:)+col3*dir;
		cntfP(faces(f,:))=cntfP(faces(f,:))+1;
		normsF(f,:)=dir;
		cntfF(f)=cntfF(f)+1;
  end
  [normsP(:,:),surfP]=normalize_rows_3D(normsP);
	[normsF(:,:),surfF]=normalize_rows_3D(normsF);

  %surfP are the surfaces areas around each point
	surfP(:)=surfP./cntfP;
  %surfF are the surface areas of each face
	surfF(:)=surfF./cntfF;
  % well we basically get them for free :)
end

%  Here we do differential geometries
% We compute two curvatures at each point
% Cp is the curvature on the short axis
% Cs is the curvature on the long axis (along the cell length)
function [CP,CS]=compute_curvatures(points,faces,normsP,np,nf)
  % initiatlization
  poins(np).lCs=[];
  gradN=zeros(3,3);
  NP=zeros(3,1);
  NS=zeros(3,1);
  dirx=ones(np,1)*[1 0 0];

  % we create the basis
  n_psi=normalize_rows_3D(cross(normsP,dirx,2));
  n_sss=cross(n_psi,normsP,2);

  % We need to make a list of the local curvatures (shared with neighbours)
  % and their weight (the smaller the distance, the less precise we are)
  % because we'll have to do some dirty tricks to clean up from numerical errors
  for i=1:np
  	poins(i).lCs=[];
  	poins(i).lCp=[];
  	poins(i).wCs=[];
  	poins(i).wCp=[];
  end

  % for all faces, we will compute a bunch of curvatures !
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

  % The cleanup is here :)
  [CP,CS]=get_curvs(poins,np);
end

% We average the curvatures, while getting read of extreme values
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

% Using simple geometries we can quickly compute the volume and surface of a discretized space
function [surf,vol]=get_surface_volume(points,face)
  nf=size(face,1);
  surf=0;
  vol=0;
  % Marginal gain
  A=zeros(1,3);
  B=zeros(1,3);
  C=zeros(1,3);
  for f=1:nf
      A(:)=points(face(f,1),:);
      B(:)=points(face(f,2),:);
      C(:)=points(face(f,3),:);
      surf=surf+norm(cross((B-A),(C-A)));
      vol=vol+abs(dot(A,cross(B,C)));
  end
  surf=surf/2;
  vol=vol/6.0;
end

% Computing the mean distance between connected points
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

% we chop the poor cells in thin slices and we measure their shapes
function [res]=make_analyze_slices(points,options)
  % We prepare the slicing
  %h=options.thickness_central/2.0;
  h=options.slice.thickness/2.0;
  mp=min(points(:,1));
  Mp=max(points(:,1));
  N=floor((Mp-mp-2*h)/(2*h));

  % we prepare the results
  res.pers=zeros(N,1);
  res.ares=zeros(N,1);
  res.pos=zeros(N,1);
  res.dist=zeros(N,1);
  res.latYZ=zeros(N,2);
  res.per(N).points=[0 0 0];
  xm=mp+2*h;

  % for all slices
  for i=1:N
      pos=xm+(i*h*2)-h;
      res.pos(i)=pos;
      res.dist(i)=min(abs(mp-pos),abs(Mp-pos));

      gl=logical((points(:,1)<(pos+h)).*(points(:,1)>(pos-h)));
      slice=points(gl,:);
      % computing eigenvalues
      [~,~,res.latYZ(i,:)]=princom(slice(:,2:3));
      % computing a smoother slice
      per= smooth_periodic( slice(:,2:3),options.central.spline);

      % from that we get our sweet sweet results
      per=[per,per(:,1)];
      [perim,ares,excl,circ]=compute_slice_shape(per);
      %res.circ(i)=1-sqrt(abs(res.latYZ(i,1)-res.latYZ(i,2))/sum(res.latYZ(i,:)));
      res.circ(i)=circ;
      res.ares(i)=ares;
      res.pers(i)=perim;
      res.excel(i)=excl;
      res.per(i).points=[res.pos(i)*ones(size(per,2),1) per(1:2,:)'];
  end
  res.nslice=N;

end

function disp_if_verbose(msg,options)
  if options.verbose>0
    disp(msg)
  end
end

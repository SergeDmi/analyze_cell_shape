function [ske,pts,arclen]=get_rolling_backbone(points,options)
% Computes the backbone of a discretized shell using rolling average

% If we really need, we can do it iteratively,
% because sometimes x is a bad approximation for the arclength
if isfield(options,'n_refine')
  n_refine=options.n_refine;
else
  n_refine=1;
end

% First we create segments
n_pts=size(points,1);
[xmin]=min(points(:,1));
[xmax]=max(points(:,1));
dx=options.dx;
X=(xmin+dx/2.0):dx:(xmax-dx/2.0);
nx=numel(X);
ske=zeros(nx,3);
ske(:,1)=X;

% Now we (iteratively) refine the backbone
for t=1:n_refine
  [ske,dirs]=refine_backbone(ske,points,dx);
end

[pts,arclen]=get_polifit_backbone(ske);

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

% Interpolates a curve regularly along its arclength
function ske=interpolate_skeleton(ske)
  s=get_arclength(ske);
  S=linspace(0,max(s),numel(s));
  if sum(s-S)~=0
	  x=spline(s,ske(:,1),S);
	  y=spline(s,ske(:,2),S);
	  z=spline(s,ske(:,3),S);
	  ske=[x',y',z'];
  end
end

% compute oridentation at each point of a discretized curve.
% extremities use a prior or posterior derivative
% other points use a central derivative
function [dirs]=get_directions(ske)
  sk=size(ske);
  nk=sk(1);
  dirf=diff(ske);
  dirs=zeros(sk);
  if nk>2
    dirs(2:(nk-1),:)=(dirf(1:(nk-2),:)+dirf(2:(nk-1),:))/2.0;
    dirs(1,:)=dirf(1,:);
    dirs(end,:)=dirf(end,:);
  else
    dirs(1,:)=dirf;
    dirs(2,:)=dirf;
  end
  dirs=normalize_rows_3D(dirs);
end

% Takes a backbone ske, and makes it a closer backbone of the points points
function [new_ske,new_dirs]=refine_backbone(ske,points,dx)
sk=size(ske);
sp=size(points);
col1=ones(sp(1),1);
row1=ones(1,3);
% First we interpolate to have points regularly spaced along the arclength, without using the information on the points
ske(:,:)=interpolate_skeleton(ske);
dirs=get_directions(ske);
% Now the skeleton is the rolling average of the points, around the previous backbone points
for i=1:sk(1)
  ptc=points-col1*ske(i,:);
  dists=abs(dot(ptc,col1*dirs(i,:),2));
  weights=exp(-(dists.^2)./(2.0*dx*dx));
  new_ske(i,:)=sum(points.*(weights*row1),1 )./sum( weights,1);
end
new_dirs=get_directions(new_ske);
end


function [pts,arclen]=get_polifit_backbone(ske)
% we compute the arclength along the curve

ske=ske(3:(end-2),:);

arclen=get_arclength(ske)';
% we align the curve
%[~,ske]=princom(ske);
% we compute the curvature of y,z as a function of the arclength

px=polyfit(arclen,ske(:,1),3);
py=polyfit(arclen,ske(:,2),3);
pz=polyfit(arclen,ske(:,3),3);

X=(arclen.^3)*px(1)+(arclen.^2)*px(2)+arclen*px(3)+px(4);
Y=(arclen.^3)*py(1)+(arclen.^2)*py(2)+arclen*py(3)+py(4);
Z=(arclen.^3)*pz(1)+(arclen.^2)*pz(2)+arclen*pz(3)+pz(4);

pts=[X Y Z];


end

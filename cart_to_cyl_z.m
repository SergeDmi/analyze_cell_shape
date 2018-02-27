function [ cylpts ] = cart_to_cyl_z( pts)
%Convert points in cartesian coordiates to points in cylindrical coords
% Cylinder axis is the z axis 
% m is the reference point in pts (angle=0)
%
%pts is of the format
%   x1 x2 x3 ... xn
%   y1 y2 y3 ... yn
%   z1 z2 z3 ... zn
%
% S.D. EMBL 2015

n=size(pts,2);
Vc=[0;0;1]*ones(1,n);
PtZ=pts(3,:); 
DPz=[pts(1:2,:);zeros(1,n)];
% projection on the plane
% Computing the angle, starting from the reference point
DOA=[1;0;0]*ones(1,n);   % Ref point
DPr=sqrt(sum(DPz.^2,1));   % Radial distance 
DV=cross(DOA,DPz,1);       % Cross product with DOA
Th=real(acos((dot(DOA,DPz,1)./(sqrt(sum(DOA.^2,1)).*DPr)))); % polar angle
% Now attributing a sign to the polar angle using the plane orientation
neg=logical(dot(DV,Vc,1)<0);
Th(neg)=-Th(neg); 

cylpts=[Th;DPr;PtZ];

end


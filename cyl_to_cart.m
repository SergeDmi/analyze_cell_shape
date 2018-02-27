function [ pts ] = cyl_to_cart( pts)
%Convert points in cylindrical coordiates to points in cartesian coords
%
%   pts are organized as :
% t1 t2 t3 ... tn
% r1 r2 r3 ... rn
% z1 z2 z3 ... zn
%
%Where  t is the angle in plane
%       r is the distance in plane
%       z is the z (out of plane_ness)
%
% S.D. EMBL 2015

pts=[pts(2,:).*cos(pts(1,:));pts(2,:).*sin(pts(1,:));pts(3,:)];

end


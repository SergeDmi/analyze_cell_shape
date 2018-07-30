function [  ] = plot_pombe( pombe,np,dx )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
pts=pombe.points;
if nargin<2
	np=size(pts,1);
end
if nargin<3
	dx=5;
end
norms=pombe.normals;

if numel(np)==1
	bli=1:np;
else
	bli=np;
end
figure
hold all

scatter3(pts(bli,1),pts(bli,2),pts(bli,3))
size([pts(bli,1) pts(bli,1)+dx*norms(bli,1)])
%plot3([pts(bli,1) pts(bli,1)+dx*norms(bli,1)],[pts(bli,2) pts(bli,2)+dx*norms(bli,2)],[pts(bli,3) pts(bli,3)+dx*norms(bli,3)]);
for i=1:np
	
	plot3([pts(i,1) pts(i,1)+dx*norms(i,1)],[pts(i,2) pts(i,2)+dx*norms(i,2)],[pts(i,3) pts(i,3)+dx*norms(i,3)]);
end
axis equal
end


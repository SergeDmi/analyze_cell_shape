function [shape] = import_ply(filename)
% Just a wrapper for plyread
% shape is an object with points (i.e. vertex)
%                         normals
%                         faces (triangles of point indexes)

[d,~] = plyread(filename);
if isfield(d.face,'vertex_index')
	ff = d.face.vertex_index;
elseif isfield(d.face,'vertex_indices')
	ff = d.face.vertex_indices;
end
nf = length(ff);
shape.faces=zeros(nf,3);
for i=1:nf
    % Converting 0-based count into 1-based count
    shape.faces(i,:) = ff{i}+1;
end

shape.points=[d.vertex.x,d.vertex.y,d.vertex.z];
shape.normals=[d.vertex.nx,d.vertex.ny,d.vertex.nz];

end

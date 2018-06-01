function [vec,ix]=exclude_extremes(v,r)
s=size(v);
m=s(2);
if s(1)>m
	v=v';
	m=s(1);
end

ls=std(((ones(m,m)-eye(m,m)).*(ones(m,1)*v)),0,2);
[~,ix]=sort(-ls);
ix=ix(1:(end-r));
vec=v(ix);

end
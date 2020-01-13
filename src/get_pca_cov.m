function [rotmat, rotated, eigenvalues]=get_pca_cov(data,axis)
% PCA analysis of data
%   data is a matrix
%   axis is the axis of data repetition
%  e.g. for a set of 3D points data , size(data)=3,N    , then axis=2 (default)
%  e.g. for a set of 3D points data , size(data)=N,3    , then axis=1

if nargin<2
  axis=2;
end

if axis==1
  data=data';
end

[n_rows,n_cols]=size(data);

% Centering data
center=mean(data')';
centered=data-repmat(center,1,n_cols);

% Getting covariance matrix M
M=(1 / (n_cols-1))*centered*centered';

% We get the rotation matrix from the eigenvectors...
[eigvectors,eigenvalues]=eig(M);
rotmat=eigvectors';
rotated=rotmat*data;

if axis==1
  rotated=rotated';
  rotmat=rotmat';
end

end

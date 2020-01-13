function [rotmat, rotated, eigenvalues]=get_pca_cov(data,axis)
% PCA analysis of data
%   data is a matrix
%   axis is the axis of data repetition
%  e.g. for a set of 3D points data , size(data)=3,N    , then axis=2 (default)
%  e.g. for a set of 3D points data , size(data)=N,3    , then axis=1
%
% One should definitely use MATLAB's princom
% But this way, no licence problem
%
% Serge Dmitrieff, 2020
% www.biophysics.fr

if nargin<2
  axis=2;
end

if axis==1
  data=data';
end

[~,n_cols]=size(data);

% Centering data
center=mean(data')';
centered=data-repmat(center,1,n_cols);

% Getting covariance matrix M
M=(1.0/(n_cols-1))*centered*centered';

% We get the rotation matrix from the eigenvectors...
[eigenvectors,eigenvalues]=eig(M);
rotmat=eigenvectors';
rotated=rotmat*data;

if axis==1
  rotated=rotated';
  rotmat=rotmat';
  eigenvalues=eigenvalues';
end

end

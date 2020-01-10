function [ L ] = seglength( Pts )
%Compute length of a chain of points Pts
%   Pts is a row of points
%   each point can be a vector (column)
% Ex, where Pts has 5 points in 2D
% Pts =  [ [0 1 1 2 2] ;
%          [0 0 1 0 1] ]
%
% S.D. EMBL 2015
dP=diff(Pts,1,2);
L=sum(sqrt(sum(dP.^2,1)));
end


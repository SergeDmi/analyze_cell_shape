function [ options ] = pombe_default_options()
%Default options for pombe analysis
%   
    options.verbose=1;
    options.thickness_central=10;
    
    options.spline.dt=2*pi/100;
    options.spline.npp=12;

end


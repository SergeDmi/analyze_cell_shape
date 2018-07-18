function [ options ] = pombe_default_options()
%Default options for pombe analysis
%
    options.verbose=1;

    options.centering=0;
    options.spline.dt=2*pi/100;
    options.spline.npp=8;
	  options.check_pairs=1;
    options.pixel_size=0.0714;
    options.do_slice_analysis=0;
    options.spline_dh=5*options.pixel_size;
    options.thickness_central=7*options.pixel_size;
end

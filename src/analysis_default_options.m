function [ options ] = analysis_default_options()
% Default options for shape analysis
%

  %% Options for program behaviour
  options.verbose=1;
	options.check_pairs=0;
  options.do_shape_resize=0;

  %% Options to choose what to analyze
  % these are typically slow analyses
  options.do_slice_analysis=1;
  options.do_curvature_analysis=0;
  options.do_backbone_analysis=0;

	%% Options for analysis
	options.pixel_size=1.0;
	%options.pixel_size=0.0714;
  options.centering=1;
  options.aligning=1;
	options.thickness=0.25;

  %% Options for slice analysis
  options.slice.thickness=2.0*options.thickness;
  options.slice.spline_dh=1.0*options.thickness;

	%% Options for perimeter spline
  options.central.spline.dt=2*pi/100;
  options.central.spline.npp=12;

	%% Options for backbone spline
	options.backbone.thickness=0.5;
	%options.backbone.spline.dt=options.central.spline.dt;
  %options.backbone.spline.npp=6;
	%options.backbone.end_exclusion=0.15; %exclude 8% front and back
  options.backbone.dx=options.thickness*2.0;
	options.backbone.n_refine=2;




end

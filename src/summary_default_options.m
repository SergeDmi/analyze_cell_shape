function [ options ] = summary_default_options()
%Default options for pombe analysis
%

  %% Options for program behaviour
  options.align_to_first=0;

  %% What to compare between experiments
  options.comparisons={'surface','volume','length','mean_circ'};



end

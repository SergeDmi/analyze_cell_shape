function [ results,experiments ] = summary_experiments( experiments,options)
% Compare diffent states for a set of experiments
%   e.g. : several cells are considered as several experiments
%   e.g. : several time points of the same cell are several states
%      experiments is a container for several states
%      each state is a container for points, filename, analysis
if nargin<2
    analysis_options=analysis_default_options();
    summary_options=summary_default_options();
else
  analysis_options=options.analysis_options;
  summary_options=options.summary_options;
end

% Bits of counting
Nexp=numel(experiments);
n_states=numel(experiments(1).states);

% What to compare
comparisons=summary_options.comparisons;
n_comp=numel(comparisons);
results(n_states).done=0;

for s=1:n_states
  for c=1:n_comp
    results(s).(comparisons{c}).values=zeros(Nexp,1);
  end
end




%% Looping over all experiments
for n=1:Nexp
  %% We align the shapes
	experiments(n)=align_experiments(experiments(n),options);
  checked=1;

  if summary_options.check_states
    disp('Manually checking experiments')
    disp(' press [SPACE] to keep, [q] to discard')
    checked=check_states(experiments(n));
  end


	if checked>0
		%experiments(n)=compare_shapes(experiments(n));
    for s=1:n_states
      % For each state, we analyze the shape
      experiments(n).states(s).analysis=analyze_shape(experiments(n).states(s).shape,analysis_options);

      for c=1:n_comp
        % For each state, we store what we are supposed to compare between stages
        results(s).(comparisons{c}).values(n)=experiments(n).states(s).analysis.(comparisons{c});
      end
    end
  end
	%experiments(n).pre_analysis=analyze_shape(experiments(n).pre_shape);
	%experiments(n).post_analysis=analyze_shape(experiments(n).post_states);

end



%% Wrap up :
for s=1:n_states
  for c=1:n_comp
    % For each state, we compute mean & std of
    % what we are supposed to compare between stages
    results(s).(comparisons{c}).mean=mean(results(s).(comparisons{c}).values);
    results(s).(comparisons{c}).std = std(results(s).(comparisons{c}).values);
  end
end



end

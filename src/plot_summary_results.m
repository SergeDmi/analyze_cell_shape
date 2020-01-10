function plot_summary_results(results,figure_ids,options)
  % Plotting results
  %   Results is a container for the stages, for the items to be compared :
  %   e.g. results(s).c.mean is the mean value of c at stage s
  summary_options=options.summary_options;
  comparisons=summary_options.comparisons;
  n_comp=numel(comparisons);
  n_figs=numel(figure_ids);
  n_states=numel(results);

  if isfield(options,'state_indices')
    state_indices=options.state_indices;
  else
    state_indices=1:n_states;
  end

  if isfield(options,'x_offset')
    state_indices=state_indices+options.x_offset;
  end


  means_states=zeros(n_states,1);
  stds_states =zeros(n_states,1);

  if n_comp~=n_figs
    error('Error : the number of figure ids and of comparison items should be the same')
  end

  for c=1:n_comp
    figure(figure_ids(c))
    for s=1:n_states
      means_states(s)=results(s).(comparisons{c}).mean;
      stds_states(s) =results(s).(comparisons{c}).std;
    end

    scatter(state_indices,means_states,options.color)
    for s=1:n_states
      plot([state_indices(s) state_indices(s)],[means_states(s)-stds_states(s)/2, means_states(s)+stds_states(s)/2],options.color)
    end
  end
end

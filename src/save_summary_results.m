function save_summary_results(results,options)
  % save results
  %   Results is a container for the stages, for the items to be compared :
  %   e.g. results(s).c.mean is the mean value of c at stage s
  summary_options=options.summary_options;
  comparisons=summary_options.comparisons;
  n_comp=numel(comparisons);
  n_states=numel(results);

  if isfield(options,'state_indices')
    state_indices=options.state_indices;
  else
    state_indices=1:n_states;
  end

  if isfield(options,'save_prefix')
    save_prefix=options.save_prefix;
  else
    save_prefix='results';
  end

  %% Saving !
  for c=1:n_comp
    % We need number of exps for detailed values
    n_exps=size(results(1).(comparisons{c}).values,1);
    values=zeros(n_exps,n_states);
    % this is to save average and std
    mean_std_states=zeros(n_states,2);
    for s=1:n_states
      values(:,s)=results(s).(comparisons{c}).values;
      mean_std_states(s,:)=[results(s).(comparisons{c}).mean,results(s).(comparisons{c}).std];
    end
    name=[save_prefix '_values_' comparisons{c} '.csv'];
    write_csv_with_names(name,values,[],state_indices)
    % now saving averaged values
    name=[save_prefix '_averaged_' comparisons{c} '.csv'];
    write_csv_with_names(name,mean_std_states,{'Mean','Std dev'},state_indices)


  end

end

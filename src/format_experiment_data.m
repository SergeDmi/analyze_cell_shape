function [experiments]=format_experiment_data(state_names,exp_names)
% prepares a properly formated set of experiments and states


Nexp=numel(exp_names);
n_states=numel(state_names);

experiments(Nexp).states(n_states).plyname='';

for count=1:Nexp
  for s=1:n_states
    name= [exp_names{count} state_names{s} '.ply'];
    experiments(count).states(s).name = name;
    experiments(count).states(s).shape = import_ply(name);
  end
end

end

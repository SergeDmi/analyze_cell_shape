function figure_ids=initiate_plot_from_names(names,xrange)
% Just creates a few plots 
n_names=numel(names);
figure_ids=zeros(n_names,1);

for n=1:n_names
  figure_ids(n)=figure();
  name=names{n};

  xlabel('Stage');
  ylabel(name);
  hold all;

  if nargin==2
    xlim(xrange);
  end
end

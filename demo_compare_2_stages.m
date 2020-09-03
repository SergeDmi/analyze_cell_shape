% Loading local file directory
startup
% Preparing options
options.summary_options=summary_default_options();
options.analysis_options=analysis_default_options();
% Adding names
stages={'1','2'};
exp_names={'thin_cleaned_sorbitol_9_','thin_cleaned_sorbitol_10_','thin_cleaned_sorbitol_11_','thin_cleaned_sorbitol_12_','thin_cleaned_sorbitol_13_','thin_cleaned_sorbitol_14_','thin_cleaned_sorbital_1_','thin_cleaned_sorbital_6_','thin_cleaned_sorbital_7_','thin_cleaned_sorbital_8_'};

% Formating data in proper structure
experiments=format_experiment_data(stages,exp_names);
% Comparing the stages
[results,experiments]=summary_experiments(experiments,options);
% Saving the results
save_summary_results(results,options)

% Nice summary plotting
comparisons=options.summary_options.comparisons;
figure_ids=initiate_plots_from_names(options.summary_options.comparisons,[0.5 4.5]);
plot_summary_results(results,figure_ids,options);

% Additional plotting to see detailed results
figure
scatter(results(1).volume.values,results(2).volume.values)
xlabel('Stage 1 volume')
ylabel('Stage 2 volume')
axis equal

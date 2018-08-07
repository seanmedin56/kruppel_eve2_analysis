function [traces, traces_krup] = gen_kruppel_eve2_data(elong_time,time_res, points_per_trace,...
                            num_traces, krup_traces, krup_ks, gen_ks, ...
                            rna_per_sec, fluo_per_rna,MS2_rise_time, ...
                            noise)
%This function returns traces with the features specified in the above
% elong_time: how long it takes an rna polymerase to traverse the gene
% time_res: How often the fluorescence is measured
% points_per_trace: How many measurements are taken
% num_traces: How many traces are taken
% avg_krup: matrix of time series at different ap positions for average
% kruppel
% krup_ks: assumed kruppel rate parameters
% rna_per_sec: Vector of how many polymerases are loaded per second for each state
% fluo_per_rna: How much fluorescence each rna produces when its fully
% loaded
% MS2_rise_time: How long it takes for a polymerase to fully fluoresce
% init_dist: The initial probability distribution for being in a given
% state

addpath('utilities/');

traces = cell([1,num_traces]);
traces_krup = cell([1, num_traces]);
rna_per_sec = [rna_per_sec(1), rna_per_sec(1), rna_per_sec(1), rna_per_sec(2)];
for i = 1:num_traces
    krup_trace = krup_traces{randi(length(krup_traces))};
    traces_krup{i} = krup_trace;
    traces{i} = gillespie_gen_krup(elong_time, time_res, points_per_trace, ...
                              krup_trace, krup_ks, gen_ks, rna_per_sec, ...
                              fluo_per_rna, MS2_rise_time,noise);
end
    
end


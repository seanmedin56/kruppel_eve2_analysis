function [traces, traces_krup] = gen_fatal_krup_data(elong_time,time_res, points_per_trace,...
                            num_traces, krup_traces, kon_krup, coop_coeff, gen_ks, ...
                            rna_per_sec, fluo_per_rna,MS2_rise_time, ...
                            noise)
%This function returns traces with the features specified in the above
% elong_time: how long it takes an rna polymerase to traverse the gene
% time_res: How often the fluorescence is measured
% points_per_trace: How many measurements are taken
% num_traces: How many traces are taken
% avg_krup: matrix of time series at different ap positions for average
% kruppel
% kon_krup: assumed kruppel on rate
% coop_coeff: number of kruppel needed to bind at once
% rna_per_sec: Vector of how many polymerases are loaded per second for each state
% fluo_per_rna: How much fluorescence each rna produces when its fully
% loaded
% MS2_rise_time: How long it takes for a polymerase to fully fluoresce
% init_dist: The initial probability distribution for being in a given
% state

addpath('utilities/');

traces = cell([1,num_traces]);
traces_krup = cell([1, num_traces]);
for i = 1:num_traces
    krup_trace = krup_traces{randi(length(krup_traces))};
    traces_krup{i} = krup_trace;
    traces{i} = gillespie_fatal_krup(elong_time, time_res, points_per_trace, ...
                              krup_trace, kon_krup, coop_coeff, gen_ks, rna_per_sec, ...
                              fluo_per_rna, MS2_rise_time,noise);
end
    
end


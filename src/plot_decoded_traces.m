addpath('utilities/');

%imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

sets_to_use = [1, 6, 8, 10, 11];
%sets_to_use = [2:5, 9, 12];
min_ap = 1;

sets_struct = trace_struct_final(ismember([trace_struct_final.setID], sets_to_use));
sets_struct = sets_struct([sets_struct.MeanAP] > min_ap);

num_bins = 30;
all_krup = [sets_struct.protein_interp];
min_krup = min(all_krup(all_krup > 0));
max_krup = max(all_krup);
krup_jump = (max_krup - min_krup) / num_bins;
krup_bins = min_krup:krup_jump:max_krup;
set_offs = cell(1, length(sets_to_use));
all_offs = cell(1, length(krup_bins));
elong_time = 8;
rise_time = 2;
thresh = 6000;

% variables I'm saving for plots
durations = [];
krup_levels = [];
time_left = [];
aps = [];

for i = 1:length(sets_to_use)
    set = sets_to_use(i);
    set_struct = sets_struct([sets_struct.setID] == set);
    set_offs{i} = cell(1, length(krup_bins));
    for trace = set_struct
        [fitted_trace, x] = fit_bursts(trace.fluo_interp, elong_time, rise_time);
        points_on = find(x > thresh);
        if length(points_on) < 2
            continue
        end
        for j = 1:(length(points_on) - 1)
            if points_on(j + 1) - points_on(j) == 1
                continue
            end
            dur = points_on(j + 1) - points_on(j) - 1;
            if dur < 0
                continue
            end
            krup_val = mean(trace.protein_interp((points_on(j) + 1):(points_on(j+1) - 1)));
            re_time = points_on(end) - points_on(j+1);
            durations = [durations dur];
            krup_levels = [krup_levels krup_val];
            time_left = [time_left re_time];
            aps = [aps trace.MeanAP];
        end
    end
end

% makes plots

figure();
scatter(krup_levels, durations, [], aps);
xlabel('krup levels');
ylabel('off duration');
title('comparing off durations and avg krup level');

figure();
scatter(krup_levels, time_left, [], aps);
xlabel('krup levels');
ylabel('time left');
title('comparing remaining time after off and avg krup level');

figure();
scatter(durations, time_left, [], aps);
xlabel('off duration');
ylabel('time left');
title('comparing remaining time after off and off durations');

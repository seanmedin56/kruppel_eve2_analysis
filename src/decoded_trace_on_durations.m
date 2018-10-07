% Looks at percent off by kruppel bin
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
set_ons = cell(1, length(sets_to_use));
all_ons = cell(1, length(krup_bins));
elong_time = 6;
rise_time = 2;
thresh = 8000;

for i = 1:length(sets_to_use)
    set = sets_to_use(i);
    set_struct = sets_struct([sets_struct.setID] == set);
    set_ons{i} = cell(1, length(krup_bins));
    for trace = set_struct
        [fitted_trace, x] = fit_bursts(trace.fluo_interp, elong_time, rise_time);
        points_off = find(x <= thresh);
        if length(points_off) < 2
            continue
        end
        for j = 1:(length(points_off) - 1)
            if points_off(j + 1) - points_off(j) == 1
                continue
            end
            dur = points_off(j + 1) - points_off(j) - 1;
            krup_val = sum(trace.protein_interp((points_off(j) + 1):(points_off(j+1) - 1)));
            krup_val = krup_val / dur;
            krup_idx = floor((krup_val - min_krup) / krup_jump) + 1;
            set_ons{i}{krup_idx} = [set_ons{i}{krup_idx} dur];
            all_ons{krup_idx} = [all_ons{krup_idx} dur];

        end
    end
end

% plots results

for i = 1:length(sets_to_use)
    mean_dur = zeros(1, length(krup_bins));
    std_dur = zeros(1, length(krup_bins));
    for j = 1:length(krup_bins)
        if isempty(set_ons{i}{j})
            mean_dur(j) = nan;
            std_dur(j) = nan;
        else
            mean_dur(j) = mean(set_ons{i}{j});
            std_dur(j) = std(set_ons{i}{j});
        end
    end
    %percent_offs(isnan(percent_offs)) = 0;
    figure();
    errorbar(krup_bins, mean_dur, std_dur);
    xlabel('kruppel');
    ylabel('time steps');
    title(['mean on duration vs krup set ' num2str(sets_to_use(i))]);
end

mean_dur = zeros(1, length(krup_bins));
std_dur = zeros(1, length(krup_bins));
for j = 1:length(krup_bins)
    if isempty(all_ons{j})
        mean_dur(j) = nan;
        std_dur(j) = nan;
    else
        mean_dur(j) = mean(all_ons{j});
        std_dur(j) = std(all_ons{j});
    end
end

figure();
errorbar(krup_bins, mean_dur, std_dur);
xlabel('kruppel');
ylabel('time steps');
title('duration on by krup level all sets');
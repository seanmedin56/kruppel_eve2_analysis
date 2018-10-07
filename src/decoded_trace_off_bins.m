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
set_counts = cell(1, length(sets_to_use));
set_offs = cell(1, length(sets_to_use));
all_counts = zeros(1, length(krup_bins));
all_offs = zeros(1, length(krup_bins));
elong_time = 8;
rise_time = 2;
thresh = 6000;

for i = 1:length(sets_to_use)
    set = sets_to_use(i);
    set_struct = sets_struct([sets_struct.setID] == set);
    set_counts{i} = zeros(1, length(krup_bins));
    set_offs{i} = zeros(1, length(krup_bins));
    for trace = set_struct
        [fitted_trace, x] = fit_bursts(trace.fluo_interp, elong_time, rise_time);
        points_on = find(x > thresh);
        if isempty(points_on)
            continue
        end
        first_point = points_on(1);
        last_point = points_on(end);
        for j = first_point:last_point
            krup_val = trace.protein_interp(j);
            krup_idx = floor((krup_val - min_krup) / krup_jump) + 1;
            set_counts{i}(krup_idx) = set_counts{i}(krup_idx) + 1;
            all_counts(krup_idx) = all_counts(krup_idx) + 1;
            if x(j) <= thresh
                set_offs{i}(krup_idx) = set_offs{i}(krup_idx) + 1;
                all_offs(krup_idx) = all_offs(krup_idx) + 1;
            end
        end
    end
end

% plots results

for i = 1:length(sets_to_use)
    percent_offs = set_offs{i} ./ set_counts{i};
    %percent_offs(isnan(percent_offs)) = 0;
    figure();
    plot(krup_bins, percent_offs);
    xlabel('kruppel');
    ylabel('percent off');
    title(['percent off by krup level set ' num2str(sets_to_use(i))]);
end

percent_offs = all_offs ./ all_counts;
%percent_offs(isnan(percent_offs)) = 0;
figure();
plot(krup_bins, percent_offs);
xlabel('kruppel');
ylabel('percent off');
title('percent off by krup level all sets');

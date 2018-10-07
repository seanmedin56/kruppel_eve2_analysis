% looks at burst magnitude as function of kruppel concentration

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
set_mags = cell(1, length(sets_to_use));
all_mags = cell(1, length(krup_bins));
elong_time = 8;
rise_time = 2;
thresh = 6000;

for i = 1:length(sets_to_use)
    set = sets_to_use(i);
    set_struct = sets_struct([sets_struct.setID] == set);
    set_mags{i} = cell(1, length(krup_bins));
    for trace = set_struct
        [fitted_trace, x] = fit_bursts(trace.fluo_interp, elong_time, rise_time);
        points_on = find(x > thresh);
        if length(points_on) < 2
            continue
        end
        for j = 2:(length(points_on) - 1)
            if points_on(j + 1) - points_on(j - 1) == 2
                continue
            end

            krup_val = trace.protein_interp(points_on(j));
            if krup_val == 0
                continue
            end
            krup_idx = floor((krup_val - min_krup) / krup_jump) + 1;
            set_mags{i}{krup_idx} = [set_mags{i}{krup_idx} x(points_on(j))];
            all_mags{krup_idx} = [all_mags{krup_idx} x(points_on(j))];

        end
    end
end

for i = 1:length(sets_to_use)
    mean_dur = zeros(1, length(krup_bins));
    std_dur = zeros(1, length(krup_bins));
    for j = 1:length(krup_bins)
        if isempty(set_mags{i}{j})
            mean_dur(j) = nan;
            std_dur(j) = nan;
        else

            mean_dur(j) = mean(set_mags{i}{j});
            std_dur(j) = std(set_mags{i}{j});
        end
    end
    %percent_offs(isnan(percent_offs)) = 0;
    figure();
    errorbar(krup_bins, mean_dur, std_dur);
    xlabel('kruppel');
    ylabel('mean eve2');
    title(['kruppel vs mean eve2 set ' num2str(sets_to_use(i))]);
end

mean_dur = zeros(1, length(krup_bins));
std_dur = zeros(1, length(krup_bins));
for j = 1:length(krup_bins)
    if isempty(all_mags{j})
        mean_dur(j) = nan;
        std_dur(j) = nan;
    else

        mean_dur(j) = mean(all_mags{j});
        std_dur(j) = std(all_mags{j});
    end
end

figure();
errorbar(krup_bins, mean_dur, std_dur);
xlabel('kruppel');
ylabel('mean eve2');
title('kruppel vs mean eve2 all sets');

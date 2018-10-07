% averages over every k time points of eve2, kruppel (prior to turn off)

% imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);
sets_to_use = [1, 6, 8, 10, 11];
set_struct = trace_struct_final(ismember([trace_struct_final.setID], sets_to_use));

k = 8;
thresh = 50000;
bins = 20;
elong_bounds = 8;

% finds highest relevant kruppel values
high_krup = 0;
for i = 1:length(set_struct)
    on_idxes = find(set_struct(i).fluo_interp > thresh);
    if ~isempty(on_idxes)
        high_krup = max([set_struct(i).protein_interp(1:on_idxes(end)) high_krup]);
    end
end
krup_jump = high_krup / bins;
krup_bins = 0:krup_jump:high_krup;

eve2_vals = cell(1, length(krup_bins));

for i = 1:length(set_struct)
    tr_eve2 = set_struct(i).fluo_interp;
    on_idxes = find(tr_eve2 > thresh);
    if length(on_idxes) - 2 * elong_bounds < k
        continue
    end
    last_idx = on_idxes(end) - elong_bounds;
    first_idx = on_idxes(1) + elong_bounds;

    tr_kr = set_struct(i).protein_interp;
    for j = first_idx:k:last_idx
        eve2_val = sum(tr_eve2(j:(j+k-1))) / k;
        kr_val = sum(tr_kr(j:(j+k-1))) / k;
        kr_idx = floor(kr_val / krup_jump) + 1;
        eve2_vals{kr_idx} = [eve2_vals{kr_idx} eve2_val];
    end   
end

mean_eve2s = zeros(1, length(krup_bins));
std_eve2s = zeros(1, length(krup_bins));
for i = 1:length(eve2_vals)
    if ~isempty(eve2_vals{i})
        mean_eve2s(i) = mean(eve2_vals{i});
        std_eve2s(i) = std(eve2_vals{i});
    end
end

figure();
errorbar(krup_bins, mean_eve2s, std_eve2s);
title([num2str(k) / 3 ' minute average eve2 by krup level']);
xlabel('kruppel level');
ylabel('eve2 level');


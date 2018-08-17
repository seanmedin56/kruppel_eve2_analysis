% looks at correlations between eve2 levels and kruppel

%imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

sets_to_use = [1 4:6 8:12];

% makes scatterplot of eve2 vs kr
eve2_vals = [];
kr_vals = [];
set_struct = trace_struct_final(ismember([trace_struct_final.setID], sets_to_use));
set_struct = set_struct([set_struct.MeanAP] > 0);
for trace = set_struct
    eve2_vals = [eve2_vals trace.fluo_interp];
    kr_vals = [kr_vals trace.protein_interp];
end
figure();
scatter(kr_vals, eve2_vals);
title('Eve2 vs Kruppel values');
xlabel('kruppel levels');
ylabel('eve2 levels');
display(corrcoef(kr_vals, eve2_vals));

% makes plot of average eve2 for each kruppel bin
max_krup = max([set_struct.protein_interp]);
min_krup = min([set_struct.protein_interp]);
num_bins = 50;
krup_jump = (max_krup - min_krup) / 50;
krup_bins = min_krup:krup_jump:max_krup;
eve2_vals = cell(1, length(krup_bins));
eve2_all = cell(1, length(krup_bins));
thresh = 50000;

for trace = set_struct
    tr_eve2 = trace.fluo_interp;
    tr_kr = trace.protein_interp;
    valid_idxes = find(tr_eve2 >= thresh);
    for i = 1:length(tr_eve2)
        kr_idx = floor((tr_kr(i) - min_krup) / krup_jump) + 1;
        eve2_all{kr_idx} = [eve2_all{kr_idx} tr_eve2(i)];
    end
    if isempty(valid_idxes)
        continue
    end
    last_idx = valid_idxes(end);
    for i = 1:last_idx
        kr_idx = floor((tr_kr(i) - min_krup) / krup_jump) + 1;
        eve2_vals{kr_idx} = [eve2_vals{kr_idx} tr_eve2(i)];
    end
end

eve2_means = zeros(1, length(eve2_vals));
eve2_means_nonzero = zeros(1, length(eve2_vals));
eve2_all_means = zeros(1, length(eve2_vals));
eve2_stds = zeros(1, length(eve2_vals));
eve2_stds_nonzero = zeros(1, length(eve2_vals));
eve2_all_stds = zeros(1, length(eve2_vals));
eve2_set_stds = zeros(1, length(eve2_vals));
for i = 1:length(eve2_vals)
    if ~isempty(eve2_vals{i})
        eve2_means(i) = mean(eve2_vals{i});
        eve2_stds(i) = std(eve2_vals{i});
        eve2_nonzero = eve2_vals{i}(eve2_vals{i} > 0);
        if ~isempty(eve2_nonzero)
            eve2_means_nonzero(i) = mean(eve2_nonzero);
            eve2_stds_nonzero(i) = std(eve2_nonzero);
        end
    end
    
    if ~isempty(eve2_all{i})
        eve2_all_means(i) = mean(eve2_all{i});
        eve2_all_stds(i) = std(eve2_all{i});
    end

    
end
figure();
errorbar(krup_bins, eve2_means, eve2_stds);
title('Eve2 vs kruppel content with 0s');
xlabel('krup content');
ylabel('mean eve2 content');

figure();
errorbar(krup_bins, eve2_means_nonzero, eve2_stds_nonzero);
title('Eve2 vs kruppel content without 0s');
xlabel('krup content');
ylabel('mean eve2 content');

figure();
errorbar(krup_bins, eve2_all_means, eve2_all_stds);
title('Eve2 vs kruppel everywhere');
xlabel('krup content');
ylabel('mean eve2 content');

% compare data across sets
figure();
for i = 1:length(sets_to_use)
    specific_set = set_struct([set_struct.setID] == sets_to_use(i));
    eve2_set_vals = cell(1, length(eve2_vals));
    eve2_set_means = zeros(1, length(eve2_set_means));
    for trace = specific_set
        tr_eve2 = trace.fluo_interp;
        tr_kr = trace.protein_interp;
        valid_idxes = find(tr_eve2 >= thresh);
        if isempty(valid_idxes)
            continue
        end
        last_idx = valid_idxes(end);
        for j = 1:last_idx
            kr_idx = floor((tr_kr(j) - min_krup) / krup_jump) + 1;
            eve2_set_vals{kr_idx} = [eve2_set_vals{kr_idx} tr_eve2(j)];
        end
    end
    for j = 1:length(eve2_set_means)
        if ~isempty(eve2_set_vals{j})
            eve2_set_means(j) = mean(eve2_set_vals{j}(eve2_set_vals{j} > thresh));
        end
    end
    plot(eve2_set_means);
    hold on
end
title('set comparison');
xlabel('krup content');
ylabel('average eve2 content');



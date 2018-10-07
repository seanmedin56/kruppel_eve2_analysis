% finds how on/off duration changes as function of kruppel content

% imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

min_ap = -0.5;
trace_struct_final = trace_struct_final([trace_struct_final.MeanAP] > min_ap);
sets_to_use = [1 4:6 8:12];

half_elong = 3;
thresh = 50000;
bin_to_show = 5;
acceptable_range = 0.5;
all_krup = [trace_struct_final.protein_interp];
all_krup = all_krup(all_krup > 0);
min_krup = min(all_krup);
max_krup = max(all_krup);
num_bins = 20;
krup_step = (max_krup - min_krup) / num_bins;
krup_bins = min_krup:krup_step:(max_krup - 1);

all_on_krup = cell(1, length(krup_bins));
all_off_krup = cell(1, length(krup_bins));

for set = sets_to_use
    on_krup = cell(1, length(krup_bins));
    off_krup = cell(1, length(krup_bins));
    set_struct = trace_struct_final([trace_struct_final.setID] == set);
    for trace = set_struct
        comp_range = max(trace.protein_interp) - min(trace.protein_interp);
        comp_mean = mean(trace.protein_interp);
        if comp_range > comp_mean * acceptable_range
            continue;
        end
        below_thresh = find(trace.fluo_interp < thresh);
        for i = 1:(length(below_thresh) - 1)
            if below_thresh(i+1) - below_thresh(i) > 1
                krup_val = mean(trace.protein_interp(...
                    below_thresh(i):below_thresh(i+1)));
                krup_idx = floor((krup_val - min_krup) / krup_step) + 1;
                krup_idx = max(min(length(krup_bins), krup_idx),1);
                on_krup{krup_idx} = [on_krup{krup_idx} ...
                    below_thresh(i+1) - below_thresh(i) - 1];
            end
        end
        above_thresh = find(trace.fluo_interp > thresh);
        for i = 1:(length(above_thresh) - 1)
            if above_thresh(i+1) - above_thresh(i) > 1
                krup_val = mean(trace.protein_interp(...
                    above_thresh(i):above_thresh(i+1)));
                krup_idx = floor((krup_val - min_krup) / krup_step) + 1;
                krup_idx = max(min(length(krup_bins), krup_idx),1);
                off_krup{krup_idx} = [off_krup{krup_idx} ...
                    above_thresh(i+1) - above_thresh(i) - 1];
            end
        end
    end
    
    % adds this sets data points to the aggregate set
    for i = 1:length(krup_bins)
        all_on_krup{i} = [all_on_krup{i} on_krup{i}];
        all_off_krup{i} = [all_off_krup{i} off_krup{i}];
    end
    
    % calculates and plots set data
    on_krup_means = zeros(1, length(on_krup));
    on_krup_stds = zeros(1, length(on_krup));
    off_krup_means = zeros(1, length(off_krup));
    off_krup_stds = zeros(1, length(off_krup));
    for i = 1:length(krup_bins)
        if ~isempty(on_krup{i})
            on_krup_means(i) = mean(on_krup{i});
            on_krup_stds(i) = std(on_krup{i});
        end
        if ~isempty(off_krup{i})
            off_krup_means(i) = mean(off_krup{i});
            off_krup_stds(i) = std(off_krup{i});
        end
    end
    figure();
    errorbar(krup_bins, on_krup_means, on_krup_stds);
    hold on
    errorbar(krup_bins, off_krup_means, off_krup_stds);
    title(['on/off durations vs kruppel a.u. set: ' num2str(set)]);
    xlabel('kruppel bin');
    ylabel('number time points');
    legend('mean on time', 'mean off time');
    
    % plots distributions at specific kruppel bin
    figure();
    histogram(on_krup{bin_to_show});
    title(['on duration distribution for specific value set ' num2str(set)]);
    
    figure();
    histogram(off_krup{bin_to_show});
    title(['off duration distribution for specific value set ' num2str(set)]);
    
end

% calculates and plots aggregate data sets
on_krup_means = zeros(1, length(on_krup));
on_krup_stds = zeros(1, length(on_krup));
off_krup_means = zeros(1, length(off_krup));
off_krup_stds = zeros(1, length(off_krup));
for i = 1:length(krup_bins)
    if ~isempty(all_on_krup{i})
        on_krup_means(i) = mean(all_on_krup{i});
        on_krup_stds(i) = std(all_on_krup{i});
    end
    if ~isempty(all_off_krup{i})
        off_krup_means(i) = mean(all_off_krup{i});
        off_krup_stds(i) = std(all_off_krup{i});
    end
end
figure();
errorbar(krup_bins, on_krup_means, on_krup_stds);
hold on
errorbar(krup_bins, off_krup_means, off_krup_stds);
title('on/off durations vs kruppel a.u. all sets');
xlabel('kruppel bin');
ylabel('number time points');
legend('mean on time', 'mean off time');

% plots distributions at specific kruppel bin
figure();
histogram(all_on_krup{bin_to_show});
title('on duration distribution for specific value all sets');

figure();
histogram(all_off_krup{bin_to_show});
title('off duration distribution for specific value all sets');

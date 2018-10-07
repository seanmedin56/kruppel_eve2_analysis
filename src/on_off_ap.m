% finds how on/off duration changes as function of ap


% imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

sets_to_use = [1 4:6 8:12];

half_elong = 3;
thresh = 50000;
bin_to_show = 12;
all_ap = [trace_struct_final.MeanAP];
min_ap = min(all_ap);
max_ap = max(all_ap);
num_bins = 25;
ap_step = (max_ap - min_ap) / num_bins;
ap_bins = min_ap:ap_step:(max_ap - 1);

all_on_ap = cell(1, length(ap_bins));
all_off_ap = cell(1, length(ap_bins));

for set = sets_to_use
    on_ap = cell(1, length(ap_bins));
    off_ap = cell(1, length(ap_bins));
    set_struct = trace_struct_final([trace_struct_final.setID] == set);
    for trace = set_struct
        below_thresh = find(trace.fluo_interp < thresh);
        for i = 1:(length(below_thresh) - 1)
            if below_thresh(i+1) - below_thresh(i) > 1
                ap_val = trace.MeanAP;
                ap_idx = floor((ap_val - min_ap) / ap_step) + 1;
                ap_idx = max(min(length(ap_bins), ap_idx),1);
                on_ap{ap_idx} = [on_ap{ap_idx} ...
                    below_thresh(i+1) - below_thresh(i) - 1];
            end
        end
        above_thresh = find(trace.fluo_interp > thresh);
        for i = 1:(length(above_thresh) - 1)
            if above_thresh(i+1) - above_thresh(i) > 1
                ap_val = trace.MeanAP;
                ap_idx = floor((ap_val - min_ap) / ap_step) + 1;
                ap_idx = max(min(length(ap_bins), ap_idx),1);
                off_ap{ap_idx} = [off_ap{ap_idx} ...
                    above_thresh(i+1) - above_thresh(i) - 1];
            end
        end
    end
    
    % adds this sets data points to the aggregate set
    for i = 1:length(ap_bins)
        all_on_ap{i} = [all_on_ap{i} on_ap{i}];
        all_off_ap{i} = [all_off_ap{i} off_ap{i}];
    end
    
    % calculates and plots set data
    on_ap_means = zeros(1, length(on_ap));
    on_ap_stds = zeros(1, length(on_ap));
    off_ap_means = zeros(1, length(off_ap));
    off_ap_stds = zeros(1, length(off_ap));
    for i = 1:length(ap_bins)
        if ~isempty(on_ap{i})
            on_ap_means(i) = mean(on_ap{i});
            on_ap_stds(i) = std(on_ap{i});
        end
        if ~isempty(off_ap{i})
            off_ap_means(i) = mean(off_ap{i});
            off_ap_stds(i) = std(off_ap{i});
        end
    end
    figure();
    errorbar(ap_bins, on_ap_means, on_ap_stds);
    hold on
    errorbar(ap_bins, off_ap_means, off_ap_stds);
    title(['on/off durations vs ap set: ' num2str(set)]);
    xlabel('ap');
    ylabel('number time points');
    legend('mean on time', 'mean off time');
    
    % plots distributions at specific kruppel bin
    figure();
    histogram(on_ap{bin_to_show});
    title(['on duration distribution for specific value set ' num2str(set)]);
    
    figure();
    histogram(off_ap{bin_to_show});
    title(['off duration distribution for specific value set ' num2str(set)]);
    
end

% calculates and plots aggregate data sets
on_ap_means = zeros(1, length(on_ap));
on_ap_stds = zeros(1, length(on_ap));
off_ap_means = zeros(1, length(off_ap));
off_ap_stds = zeros(1, length(off_ap));
for i = 1:length(ap_bins)
    if ~isempty(all_on_ap{i})
        on_ap_means(i) = mean(all_on_ap{i});
        on_ap_stds(i) = std(all_on_ap{i});
    end
    if ~isempty(all_off_ap{i})
        off_ap_means(i) = mean(all_off_ap{i});
        off_ap_stds(i) = std(all_off_ap{i});
    end
end
figure();
errorbar(ap_bins, on_ap_means, on_ap_stds);
hold on
errorbar(ap_bins, off_ap_means, off_ap_stds);
title('on/off durations vs ap (sec) all sets');
xlabel('ap bin');
ylabel('number time points');
legend('mean on time', 'mean off time');

% plots distributions at specific kruppel bin
figure();
histogram(all_on_ap{bin_to_show});
title('on duration distribution for specific value all sets');

figure();
histogram(all_off_ap{bin_to_show});
title('off duration distribution for specific value all sets');
% finds how on/off duration changes as a function of time

% imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

min_ap = -0.5;
trace_struct_final = trace_struct_final([trace_struct_final.MeanAP] > min_ap);
sets_to_use = [1 4:6 8:12];

half_elong = 3;
thresh = 50000;
bin_to_show = 5;
all_times = [trace_struct_final.time_interp];
min_time = 0;
max_time = max(all_times);
num_bins = 20;
time_step = (max_time - min_time) / num_bins;
time_bins = min_time:time_step:(max_time - 1);

all_on_time = cell(1, length(time_bins));
all_off_time = cell(1, length(time_bins));

for set = sets_to_use
    on_time = cell(1, length(time_bins));
    off_time = cell(1, length(time_bins));
    set_struct = trace_struct_final([trace_struct_final.setID] == set);
    for trace = set_struct
        below_thresh = find(trace.fluo_interp < thresh);
        for i = 1:(length(below_thresh) - 1)
            if below_thresh(i+1) - below_thresh(i) > 1
                time_val = trace.time_interp(below_thresh(i+1));
                time_idx = floor((time_val - min_time) / time_step) + 1;
                time_idx = max(min(length(time_bins), time_idx),1);
                on_time{time_idx} = [on_time{time_idx} ...
                    below_thresh(i+1) - below_thresh(i) - 1];
            end
        end
        above_thresh = find(trace.fluo_interp > thresh);
        for i = 1:(length(above_thresh) - 1)
            if above_thresh(i+1) - above_thresh(i) > 1
                time_val = trace.time_interp(above_thresh(i+1));
                time_idx = floor((time_val - min_time) / time_step) + 1;
                time_idx = max(min(length(time_bins), time_idx),1);
                off_time{time_idx} = [off_time{time_idx} ...
                    above_thresh(i+1) - above_thresh(i) - 1];
            end
        end
    end
    
    % adds this sets data points to the aggregate set
    for i = 1:length(time_bins)
        all_on_time{i} = [all_on_time{i} on_time{i}];
        all_off_time{i} = [all_off_time{i} off_time{i}];
    end
    
    % calculates and plots set data
    on_time_means = zeros(1, length(on_time));
    on_time_stds = zeros(1, length(on_time));
    off_time_means = zeros(1, length(off_time));
    off_time_stds = zeros(1, length(off_time));
    for i = 1:length(time_bins)
        if ~isempty(on_time{i})
            on_time_means(i) = mean(on_time{i});
            on_time_stds(i) = std(on_time{i});
        end
        if ~isempty(off_time{i})
            off_time_means(i) = mean(off_time{i});
            off_time_stds(i) = std(off_time{i});
        end
    end
    figure();
    errorbar(time_bins, on_time_means, on_time_stds);
    hold on
    errorbar(time_bins, off_time_means, off_time_stds);
    title(['on/off durations vs time (sec) set: ' num2str(set)]);
    xlabel('time');
    ylabel('number time points');
    legend('mean on time', 'mean off time');
    
    % plots distributions at specific kruppel bin
    figure();
    histogram(on_time{bin_to_show});
    title(['on duration distribution for specific value set ' num2str(set)]);
    
    figure();
    histogram(off_time{bin_to_show});
    title(['off duration distribution for specific value set ' num2str(set)]);
    
end

% calculates and plots aggregate data sets
on_time_means = zeros(1, length(on_time));
on_time_stds = zeros(1, length(on_time));
off_time_means = zeros(1, length(off_time));
off_time_stds = zeros(1, length(off_time));
for i = 1:length(time_bins)
    if ~isempty(all_on_time{i})
        on_time_means(i) = mean(all_on_time{i});
        on_time_stds(i) = std(all_on_time{i});
    end
    if ~isempty(all_off_time{i})
        off_time_means(i) = mean(all_off_time{i});
        off_time_stds(i) = std(all_off_time{i});
    end
end
figure();
errorbar(time_bins, on_time_means, on_time_stds);
hold on
errorbar(time_bins, off_time_means, off_time_stds);
title('on/off durations vs time (sec) all sets');
xlabel('time bin');
ylabel('number time points');
legend('mean on time', 'mean off time');

% plots distributions at specific kruppel bin
figure();
histogram(all_on_time{bin_to_show});
title('on duration distribution for specific value all sets');

figure();
histogram(all_off_time{bin_to_show});
title('off duration distribution for specific value all sets');

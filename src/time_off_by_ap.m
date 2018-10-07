% plots distribution of off times for certain aps for certain sets (along
% with kruppel values at each time point)

% imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

sets_to_use = [1];
aps_to_use = 3:7;
ap_width = 0.5;
time_interv = 180;
thresh = 50000;

for set = sets_to_use
    set_struct = trace_struct_final([trace_struct_final.setID] == set);
    for ap = aps_to_use
        ap_struct = set_struct(abs([set_struct.MeanAP] - ap) < ap_width);
        min_time = min([ap_struct.time_interp]);
        max_time = max([ap_struct.time_interp]);
        time_bins = min_time:time_interv:max_time;
        off_times = zeros(1, length(time_bins));
        krup_vals = cell(1, length(time_bins));
        for trace = ap_struct
            on_idxes = find(trace.fluo_interp > thresh);
            if isempty(on_idxes) || on_idxes(end) == length(trace.fluo_interp)
                continue
            end
            stop_time = trace.time_interp(on_idxes(end));
            time_idx = floor((stop_time - min_time) / time_interv) + 1;
            time_idx = min(time_idx, length(time_bins));
            off_times(time_idx) = off_times(time_idx) + 1;
            for i = 1:(length(krup_vals) - 1)
                idxes = find(trace.time_interp > time_bins(i) & ...
                    trace.time_interp < time_bins(i + 1));
                krup_vals{i} = [krup_vals{i} trace.protein_interp(idxes)];
            end
        end
        krup_means = zeros(1, length(krup_vals));
        for i = 1:(length(krup_vals) - 1)
            krup_means(i) = mean(krup_vals{i});
        end
        figure();
        yyaxis left
        plot(time_bins, off_times);
        yyaxis right
        plot(time_bins, krup_means);
        title(['off times at ap ' num2str(ap) ' of set ' num2str(set)]);
    end
end


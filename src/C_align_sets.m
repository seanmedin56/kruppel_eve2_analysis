% aligns trace's ap and time positions (and also scales the kruppel values)
% based on kruppel concentration
addpath('utilities/');
project = 'Kruppel_eve2_pass1'; %Project Identifier
dat_path = ['../dat/' project '/']; % data mat directory

%loads trace_struct_final from B
load([dat_path 'inference_traces_Kruppel_eve2_pass1_dT20.mat']);
new_final = trace_struct_final;
sets_to_use = [1:6 8:12];
alignment_scalers = ones(1, length(sets_to_use));
ap_offsets = zeros(1, length(sets_to_use));
time_offsets = zeros(1, length(sets_to_use));
set_mats = cell(1, length(sets_to_use));
fit_mats = cell(1, length(sets_to_use));
ap_bins = -1.5:0.5:3;
time_bin = 120;

% compiles each set into usable form
for s = 1:length(sets_to_use)
    set = sets_to_use(s);
    set_struct = trace_struct_final([trace_struct_final.setID] == set);
    min_time = min([set_struct.time_interp]) + 600;
    max_time = max([set_struct.time_interp]);
    time_bins = min_time:time_bin:(max_time - 1);
    mat_res = zeros(length(ap_bins), length(time_bins));
    for i = 1:length(ap_bins)
        ap = ap_bins(i);
        for j = 1:length(time_bins)
            time = time_bins(j);
            rel_traces = set_struct(abs([set_struct.MeanAP] - ap) < 0.25);
            all_times = [rel_traces.time_interp];
            all_fluo = [rel_traces.protein_interp];
            idxes = find(all_times - time >= 0 & all_times - time <= time_bin);
            mat_res(i,j) = mean(all_fluo(idxes));
        end
    end
    set_mats{s} = mat_res;
end

% aligns sets one by one (using first set as baseline)
fit_mats{1} = set_mats{1};
for i = 2:length(set_mats)
    best_scaler = 0;
    best_error = 100000000000000;
    best_time_delay = 0;
    best_ap_delay = 0;
    for ap_delay = -4:4
        for time_delay = -10:10
            [error,scaler] = get_error_and_scaler(set_mats{i}, ...
                fit_mats, time_offsets(1:(i-1)), ...
                ap_offsets(1:(i-1)), ap_delay, time_delay);

            if error < best_error
                best_error = error;
                best_scaler = scaler;
                best_time_delay = time_delay;
                best_ap_delay = ap_delay;
            end
        end
    end
    alignment_scalers(i) = best_scaler;
    fit_mats{i} = set_mats{i} * best_scaler;
    time_offsets(i) = best_time_delay;
    ap_offsets(i) = best_ap_delay;
    idxes = find([new_final.setID] == sets_to_use(i));
    for idx = idxes
        new_final(idx).time_interp = new_final(idx).time_interp + time_bin ...
            * best_time_delay;
        new_final(idx).MeanAP = new_final(idx).MeanAP + 0.5 * best_ap_delay;
        new_final(idx).protein_interp = new_final(idx).protein_interp * best_scaler;
    end
end

% readjusts times so that all times > 0
time_adjust = min([new_final.time_interp]);
for i = 1:length(new_final)
    new_final(i).time_interp = new_final(i).time_interp - time_adjust;
end
trace_struct_final = new_final;
save([dat_path 'adjust_traces.mat'], 'trace_struct_final', 'time_offsets', ...
    'ap_offsets', 'alignment_scalers');



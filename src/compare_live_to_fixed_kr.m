addpath('utilities/');
% finds ap position/time/scalar multiplier in fixed tissue experiment 
% kruppel data that best matches the live kruppel data for each set, calculates 
% the error and plots them side by side

% imports data
sets_to_use = 1:3;
live_data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
fig_path = '../fig/Kruppel_eve2_pass1/';
mkdir(fig_path);
load(live_data_path);
fixed_data_path = '../dat/mov_raw_dorsal_111211_Kni_Kr_Gt_Hb_orer_slideA.mat';
load(fixed_data_path);

ap_step = 1;
time_step = 180; % must be >= 60 seconds

% puts fixed kruppel data into usable form
KrMov_final = KrMov_final';
ap_bin_bounds = 200:10 * ap_step:600;
fixed_kr_mat = zeros(length(ap_bin_bounds), 3600 / time_step);
ap_idx = 1;
for ap = ap_bin_bounds
    start = 1;
    time_idx = 1;
    while start < 60
        fixed_kr_mat(ap_idx, time_idx) = sum(sum(KrMov_final(ap:ap + 10*ap_step-1,start:start+2)));
        start = start + 3;
        time_idx = time_idx + 1;
    end
    ap_idx = ap_idx + 1;
end

% get plots of time trends of for each ap position for every set

ap_bins = -3:ap_step:9;
time_bins = 0:time_step:3000;
consolidated_fluo = struct;
idx = 1;
for set = sets_to_use
    trace_struct_set = trace_struct_final([trace_struct_final.setID] == set);
    for ap = ap_bins
        ap_traces_set = trace_struct_set(abs([trace_struct_set.MeanAP] ...
            - ap) <= ap_step / 2);
        fluo_by_time = zeros(1, length(time_bins));
        std_fluo_by_time = zeros(1, length(time_bins));
        time_idx = 1;
        for time = time_bins
            fluos = [];
            consolidated_fluo(idx).ap = ap;
            consolidated_fluo(idx).set = set;
            for trace = ap_traces_set
                idxes = find(abs(trace.time_interp - time) <= time_step / 2);
                fluos = [fluos trace.protein_interp(idxes)];
            end
            if ~isempty(fluos)
                fluo_by_time(time_idx) = mean(fluos);
                std_fluo_by_time(time_idx) = std(fluos);
            end
            time_idx = time_idx + 1;
        end
        consolidated_fluo(idx).fluo = fluo_by_time;
        consolidated_fluo(idx).std_fluo = std_fluo_by_time;
        idx = idx + 1;
    end
    
    % finds best time/ap fit in fixed kruppel data
    
    % but first forms matrix of live kruppel data
    live_kr_mat = zeros(length(ap_bins),length(time_bins));
    ap_idx = 1;
    max_time = length(time_bins);
    min_time = 1;
    for ap = ap_bins
        trace = consolidated_fluo([consolidated_fluo.ap] == ap ...
            & [consolidated_fluo.set] == set);
        live_kr_mat(ap_idx,1:end) = trace.fluo;
        nonzero = find(trace.fluo ~= 0);
        max_time = min(nonzero(end), max_time);
        min_time = max(nonzero(1), min_time);
        ap_idx = ap_idx + 1;
    end
    live_kr_mat = live_kr_mat(1:end, min_time:max_time);
    [best_mat,best_scalar,row,col] = get_best_scaled_sub_mat(live_kr_mat, ...
        fixed_kr_mat);
    display(row);
    display(col);
    start_ap = 4;
    for ap = [-3, 0, 3, 6, 9]
        figure();
        plot(live_kr_mat(ap + start_ap,:));
        hold on
        plot(best_mat(ap + start_ap,:));
        title(['Set ' num2str(set) '. AP ' num2str(ap)]);
        xlabel('time steps');
        ylabel('adjusted fluorescence');
        legend('live', 'fixed');
    end
    
end



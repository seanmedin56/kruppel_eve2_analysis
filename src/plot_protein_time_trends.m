% plots ap/time trends of kruppel data

sets_to_use = [1 4:6 8:12];
%data_path = '../dat/Kruppel_eve2_pass1/adjust_traces.mat';
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
fig_path = '../fig/Kruppel_eve2_pass1/';
mkdir(fig_path);
load(data_path);

% makes 3D scatter plot of data if there are 5 or less data sets
if length(sets_to_use) <= 5
    figure();
    leg = cell(1, length(sets_to_use));
    idx = 1;
    for set = sets_to_use
        aps = [];
        fluos = [];
        times = [];
        trace_struct_set = trace_struct_final([trace_struct_final.setID] == set);
        for trace = trace_struct_set
            aps = [aps ones(1,length(trace.protein_interp)) * trace.MeanAP];
            times = [times trace.time_interp];
            fluos = [fluos trace.protein_interp];
        end
        scatter3(aps, times,fluos);
        leg{idx} = ['set ' num2str(set)]; 
        idx = idx + 1;
        hold on
    end
    title('Protein Fluorescence Scatter');
    legend(leg);
    xlabel('ap');
    ylabel('time');
    zlabel('fluorescence');
    savefig([fig_path, 'fluo_scatter.fig']);
end

% get plots of time trends of for each ap position for every set
ap_step = 1;
time_step = 180;
ap_bins = -5:ap_step:9;
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
            fluo_by_time(time_idx) = mean(fluos);
            std_fluo_by_time(time_idx) = std(fluos);
            time_idx = time_idx + 1;
        end
        consolidated_fluo(idx).fluo = fluo_by_time;
        consolidated_fluo(idx).std_fluo = std_fluo_by_time;
        idx = idx + 1;
    end
end

for ap = 0
    figure();
    relevant_sets = consolidated_fluo([consolidated_fluo.ap] == ap);
    leg = cell(1, length(relevant_sets));
    for i = 1:length(relevant_sets)
        fluo_by_time = relevant_sets(i).fluo;
        std_fluo_by_time = relevant_sets(i).std_fluo;
        errorbar(time_bins, fluo_by_time, std_fluo_by_time);
        leg{i} = ['set ' num2str(relevant_sets(i).set)];
        hold on
    end
    title(['Fluorescence for ap ' num2str(ap)])
    xlabel('time');
    ylabel('average fluorescence');
    legend(leg);
    savefig([fig_path 'time_trends_ap_' num2str(ap)]);
end

% plots mean and standard deviation fluo vs time among sets
for ap = -3:3:9
    figure();
    relevant_sets = consolidated_fluo([consolidated_fluo.ap] == ap);
    fluo_means = zeros(1, length(time_bins));
    fluo_stds = zeros(1, length(time_bins));
    fluo_maxes = zeros(1, length(time_bins));
    fluo_mins = zeros(1, length(time_bins));
    for i = 1:length(fluo_means)
        fluo_vals = [];
        for j = 1:length(relevant_sets)
            if relevant_sets(j).fluo(i) > 0
                fluo_vals = [fluo_vals relevant_sets(j).fluo(i)];
            end
        end
        if ~isempty(fluo_vals)
            fluo_means(i) = mean(fluo_vals);
            fluo_stds(i) = std(fluo_vals);
            fluo_maxes(i) = max(fluo_vals);
            fluo_mins(i) = min(fluo_vals);
        end
    end
    errorbar(time_bins, fluo_means, fluo_stds);
    display(['ap ' num2str(ap)]);
    rel_std = sort(fluo_stds ./ fluo_means);
    display(rel_std);
    rel_range = sort((fluo_maxes - fluo_mins) ./ fluo_means);
    display(rel_range);
    
    hold on
    plot(time_bins, fluo_maxes);
    hold on
    plot(time_bins, fluo_mins);
    title(['Fluorescence for ap ' num2str(ap)])
    xlabel('time');
    ylabel('average fluorescence');
end


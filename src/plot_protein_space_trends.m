% plots space trends at various time points

% imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

sets_to_use = [1:6 8:12];

% get plots of space trends at each specified time point
half_time_spread = 150;
time_points = 600:300:1800;
ap_points = -6:0.5:12;
half_ap_spread = 0.25;
consolidated_fluo = struct;
idx = 1;
for set = sets_to_use
    trace_struct_set = trace_struct_final([trace_struct_final.setID] == set);
    for i = 1:length(ap_points)
        ap = ap_points(i);
        ap_traces_set = trace_struct_set(abs([trace_struct_set.MeanAP] ...
            - ap) <= half_ap_spread);
        time_idx = 1;
        for j = 1:length(time_points)
            fluos = [];
            time = time_points(j);
            consolidated_fluo(idx + j - 1).time = time;
            consolidated_fluo(idx + j - 1).set = set;
            for trace = ap_traces_set
                idxes = find(abs(trace.time_interp - time) <= half_time_spread);
                fluos = [fluos trace.protein_interp(idxes)];
            end
            if ~isempty(fluos)
                fluos = [fluos 0];
            end
            consolidated_fluo(idx + j - 1).fluo(i) = mean(fluos);
            consolidated_fluo(idx + j - 1).std_fluo(i) = std(fluos);
        end
    end
    idx = idx + length(time_points);
end

% plot trends

for time = time_points
    figure();
    relevant_sets = consolidated_fluo([consolidated_fluo.time] == time);
    leg = cell(1, length(relevant_sets));
    for i = 1:length(relevant_sets)
        fluo_by_time = relevant_sets(i).fluo;
        std_fluo_by_time = relevant_sets(i).std_fluo;
        errorbar(ap_points, fluo_by_time, std_fluo_by_time);
        leg{i} = ['set ' num2str(relevant_sets(i).set)];
        hold on
    end
    title(['Fluorescence at ' num2str(time / 60) ' minutes'])
    xlabel('ap bin');
    ylabel('average fluorescence');
    legend(leg);
end

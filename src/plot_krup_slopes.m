% plots average kruppel slope over time by ap

%imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);
sets_to_use = [6];

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
                if length(idxes) > 1
                    fluos = [fluos diff(trace.protein_interp(idxes))];
                end
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

for ap = -3:3:9
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
end
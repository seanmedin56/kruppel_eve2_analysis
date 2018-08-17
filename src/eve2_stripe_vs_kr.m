% plots average eve2 intensity at center of stripe and kr concentration on
% same graph

% imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

sets_to_use = [1 4:6 8:12];

aps = -3:1:6;

for ap = aps
    mean_eve2s = [];
    mean_krs = [];
    figure();
    for set = sets_to_use
        set_struct = trace_struct_final([trace_struct_final.setID] == set);
        set_struct = set_struct(abs([set_struct.MeanAP] - ap) < 0.5);
        time_min = 900;
        time_max = 2400;
        eve2 = [];
        kr = [];
        for trace = set_struct
            idxes = find([trace.time_interp] > time_min & [trace.time_interp] ...
                < time_max);
            if ~isempty(idxes)
                eve2 = [eve2 trace.fluo_interp(idxes)];
                kr = [kr trace.protein_interp(idxes)];
            end
        end
        mean_eve2s = [mean_eve2s mean(eve2)];
        mean_krs = [mean_krs mean(kr)];
    end
    scatter(mean_eve2s, mean_krs);
    xlabel('average eve2');
    ylabel('average kr');
    title(['looking at eve2 and kr at center of stripe ap ' num2str(ap)]);
    display(corrcoef(mean_eve2s, mean_krs));
end
% compares average on eve2 values to kruppel levels among the different
% sets

% imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

sets_to_use = [2:5 9 12];

min_ap = 0;
max_ap = 9;
half_elong = 3;
thresh = 50000;

%purges traces that are higher than the minimum ap
trace_struct_final = trace_struct_final([trace_struct_final.MeanAP] > min_ap);
trace_struct_final = trace_struct_final([trace_struct_final.MeanAP] < max_ap);

all_eve2 = [];
all_kr = [];
for set = sets_to_use
    set_struct = trace_struct_final([trace_struct_final.setID] == set);
    set_eve2 = [];
    set_kr = [];
    for trace = set_struct
        below_thresh = [0 find(trace.fluo_interp < thresh) ...
            length(trace.fluo_interp) + 1];
        for i = 1:(length(below_thresh) - 1)
            if below_thresh(i+1) - below_thresh(i) > 2 * half_elong
                new_eve2 = mean(trace.fluo_interp(below_thresh(i) + ...
                    half_elong:below_thresh(i+1) - half_elong));
                new_kr = mean(trace.protein_interp(below_thresh(i) + ...
                    half_elong:below_thresh(i+1) - half_elong));
                set_eve2 = [set_eve2 new_eve2];
                set_kr = [set_kr new_kr];
            end
        end
    end
    figure();
    scatter(set_eve2, set_kr);
    xlabel('average on eve2 fluorescence');
    ylabel('average kruppel fluorescence');
    coeff_mat = corrcoef(set_eve2,set_kr);
    title(['eve2 vs kruppel set ' num2str(set) 'and corr: ' num2str(coeff_mat(1,2))]);
    
    all_eve2 = [all_eve2 set_eve2];
    all_kr = [all_kr set_kr];
end

% scatters everything
figure();
% idxes = find(all_eve2 < 200000);
% all_eve2 = all_eve2(idxes);
% all_kr = all_kr(idxes);
scatter(all_eve2, all_kr);
xlabel('average on eve2 fluorescence');
ylabel('average kruppel fluorescence');
coeff_mat = corrcoef(all_eve2,all_kr);
title(['eve2 vs kruppel all sets and corr: ' num2str(coeff_mat(1,2))]);

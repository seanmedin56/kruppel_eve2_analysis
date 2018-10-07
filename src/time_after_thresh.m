% looks at distribution of off times once threshhold value of kruppel is
% reached

%imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

sets_to_use = [10];
krup_thresh = 6000000;
eve2_thresh = 50000;

set_struct = trace_struct_final(ismember([trace_struct_final.setID], sets_to_use));
set_struct = set_struct([set_struct.MeanAP] > 1);
time_after = [];
num_exceptions = 0;

for trace = set_struct
    on_idxes = find(trace.fluo_interp > eve2_thresh);
    if isempty(on_idxes) || on_idxes(end) == length(trace.fluo_interp)
        continue
    end
    off_idx = on_idxes(end);
    past_thresh = find(trace.protein_interp >= krup_thresh);
    if isempty(past_thresh)
        num_exceptions = num_exceptions + 1;
        continue
    end
    first_past = past_thresh(1);
    time_after = [time_after off_idx - first_past];
end
figure();
histogram(time_after);
title('Time off after thresh');



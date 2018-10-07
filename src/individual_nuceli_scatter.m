% looks at average eve2 for each nuclei in center of stripe

%import data
% imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

set = 9;
ap = 3;
half_width = 3;

set_struct = trace_struct_final([trace_struct_final.setID] == set);
set_struct = set_struct(abs([set_struct.MeanAP] - ap) < half_width);

eve2_vals = [];
nonzero_eve2_vals = [];
kr_vals = [];


for trace = set_struct
    nonzero = find(trace.fluo_interp > 0);
    eve2_vals = [eve2_vals mean(trace.fluo_interp)];
    if isempty(nonzero)
        nonzero_eve2_vals = [nonzero_eve2_vals 0];
    else
        nonzero_eve2_vals = [nonzero_eve2_vals mean(trace.fluo_interp(nonzero))];
    end
    kr_vals = [kr_vals mean(trace.protein_interp)];
end
figure();
scatter(eve2_vals, kr_vals);
title('all eve2 vals');
xlabel('eve2');
ylabel('kruppel');

figure();
scatter(nonzero_eve2_vals, kr_vals);
title('nonzero eve2 vals');
xlabel('eve2');
ylabel('kruppel');

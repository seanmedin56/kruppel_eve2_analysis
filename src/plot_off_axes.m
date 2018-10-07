% compares different values when a nuclei turns off

%imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);
sets_to_use = [1,4,6,10,11];
thresh = 50000;
assumed_elong = 8;

% collects off data for each nuclei
set_struct = trace_struct_final(ismember([trace_struct_final.setID], ...
    sets_to_use));
set_struct = set_struct([set_struct.MeanAP] > 2);
off_data = [];

for trace = set_struct
    on_idxes = find(trace.fluo_interp > thresh);
    if length(on_idxes) < 2 || on_idxes(end) == length(trace.fluo_interp)
        continue;
    end
    off_tr.time = trace.time_interp(on_idxes(end));
    off_tr.last_kr = trace.protein_interp(on_idxes(end));
    if on_idxes(end) > assumed_elong
        off_tr.elong_ago_kr = trace.protein_interp(on_idxes(end) ...
            - assumed_elong);
        off_tr.elong_kr_slope = trace.protein_interp(on_idxes(end)) ...
        - trace.protein_interp(on_idxes(end) - assumed_elong);
    else
        off_tr.elong_ago_kr = NaN;
        off_tr.elong_kr_slope = NaN;
    end
    if on_idxes(end) + assumed_elong < length(trace.protein_interp)
        off_tr.elong_future_slope = trace.protein_interp(on_idxes(end) + ...
            assumed_elong) - trace.protein_interp(on_idxes(end));
        off_tr.elong_future_kr = trace.protein_interp(on_idxes(end) + ...
            assumed_elong);
    else
        off_tr.elong_future_slope = NaN;
        off_tr.elong_future_kr = NaN;
    end
    off_tr.ap = trace.MeanAP;
    off_tr.last_kr_slope = trace.protein_interp(on_idxes(end)) ...
        - trace.protein_interp(on_idxes(end) - 1);
    
    off_data = [off_data off_tr];
end

non_nans_past = ~isnan([off_data.elong_ago_kr]);
non_nans_future = ~isnan([off_data.elong_future_kr]);

% plots ap vs kruppel off
figure();
scatter([off_data.ap], [off_data.last_kr]);
xlabel('ap');
ylabel('kruppel');
title('off times ap vs kr');

% plots ap vs previous off kr
figure();
scatter([off_data(non_nans_past).ap], [off_data(non_nans_past).elong_ago_kr]);
xlabel('ap');
ylabel('kruppel');
title(['off times ap vs kr delayed ' num2str(assumed_elong) ' points']);

% plots immediate kr vs immediate kr slope
figure();
scatter([off_data.last_kr], [off_data.last_kr_slope]);
xlabel('kruppel');
ylabel('kruppel change');
title('off times kruppel vs instant kruppel change');

% plots immediate kr vs elongation time kr slope
figure();
scatter([off_data(non_nans_past).last_kr], [off_data(non_nans_past).elong_kr_slope]);
xlabel('kruppel');
ylabel('kruppel change');
title('off times kruppel vs elong time kruppel change');

% plots time vs ap
figure();
scatter([off_data.time], [off_data.ap]);
xlabel('time');
ylabel('ap');
title('time vs ap off');

% plots immediate kr vs future elongation long slope
figure();
scatter([off_data(non_nans_future).last_kr], ...
    [off_data(non_nans_future).elong_future_slope]);
xlabel('kruppel');
ylabel('kruppel change');
title('off times kruppel vs future kruppel added');


% compares width of eve2 stripe (starting from center of stripe) to kruppel
% content at boundary of stripe (at 5 minute increments)


% imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

sets_to_use = [1:6 8:12];
thresh = 0.25; % percent of nuclei needing to be on in time window
eve2_thresh = 50000;

% creates scatter plot of widths and kruppel values (different colors for
% each set)
max_widths = [];
krs_at_max_widths = [];
figure();
for set = sets_to_use
    widths = [];
    krs = [];
    set_struct = trace_struct_final([trace_struct_final.setID] == set);

    ap_incr = 0.25;
    ap_spread = 0.25;
    max_time = max([set_struct.time_interp]);
    time_incr = 180;
    
    for time = 0:time_incr:(max_time - time_incr)
        ap_val = 0;
        while true
            place_struct = set_struct(abs([set_struct.MeanAP] - ap_val) < ap_spread);
            on_nucs = 0;
            total_nucs = 0;
            kr_vals = [];
            for trace = place_struct
                idxes = find(trace.time_interp > time & ...
                    trace.time_interp < time + time_incr);
                kr_vals = [kr_vals trace.protein_interp(idxes)];
                if sum(find(trace.fluo_interp(idxes) > eve2_thresh)) > 0
                    on_nucs = on_nucs + 1;
                end
                total_nucs = total_nucs + 1;
            end
            
            if on_nucs / total_nucs < thresh
                break
            end
            ap_val = ap_val + ap_incr;
        end
        widths = [widths ap_val];
        krs = [krs mean(kr_vals)];
    end
    [max_width, max_idx] = max(widths);
    max_widths = [max_widths max_width];
    krs_at_max_widths = [krs_at_max_widths krs(max_idx)];
end
scatter(max_widths, krs_at_max_widths);



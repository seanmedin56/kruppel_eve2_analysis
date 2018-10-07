% plots width shift and kruppel level at boundary as function of time

% imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

sets_to_use = [1:6 8:12];
thresh = 0.75; % what percent of nuclei need to be off to be considered boundary
eve2_thresh = 50000;
start_ap = 0;
ap_incr = 0.25;
ap_half = 0.5;
time_interv = 180; % looks at 3 minute time windows
set_krups = cell(1, length(sets_to_use));
set_times = cell(1, length(sets_to_use));
idx = 1;
for set = sets_to_use
    set_struct = trace_struct_final([trace_struct_final.setID] == set);
    ap_vals = [];
    kruppel_vals = [];
    min_time = min([set_struct.time_interp]);
    max_time = max([set_struct.time_interp]);
    times = min_time:time_interv:max_time;
    for time = times
        ap_val = 0;
        % finds last ap that meets threshhold for being on
        while true
            ap_struct = set_struct(abs([set_struct.MeanAP] - ap_val) < ap_half);
            if isempty(ap_struct)
                break
            end
            num_on = 0;
            num_tot = 0;
            krups = [];
            for trace = ap_struct
                idxes = find(trace.time_interp >= time & trace.time_interp < ...
                    time + time_interv);
                if isempty(idxes)
                    continue
                end
                if ~isempty(find(trace.fluo_interp(idxes) > eve2_thresh, 1))
                    num_on = num_on + 1;
                end
                num_tot = num_tot + 1;
                krups = [krups trace.protein_interp(idxes)];
            end
            if num_tot == 0 || 1 - (num_on / num_tot) >= thresh
                break
            end
            ap_val = ap_val + ap_incr;
        end
        ap_vals = [ap_vals ap_val];
        kruppel_vals = [kruppel_vals mean(krups)];
    end
    figure();
    yyaxis left
    plot(times, ap_vals);
    ylabel('width of stripe');
    yyaxis right
    plot(times, kruppel_vals);
    ylabel('kruppel value at edge');
    xlabel('time');
    title(['edge vs kruppel level set ' num2str(set)])
    set_krups{idx} = kruppel_vals;
    set_times{idx} = times;
    idx = idx + 1;
end

figure();
leg = {};
for idx = 1:length(sets_to_use)
    plot(set_times{idx}, set_krups{idx});
    hold on
    leg{idx} = ['set ' num2str(sets_to_use(idx))];
end
legend(leg);
xlabel('time');
ylabel('kruppel val');
title('all sets edge kruppel values');
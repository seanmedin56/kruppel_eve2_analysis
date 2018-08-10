% script looks at instantaneous and integrated kruppel level when eve2
% turns off for the last time (eve2 is off once no more spots are detected
% for more than 4 time steps)

%imports data

sets_to_use = 1:3;
min_burst_size = 5;
min_nonzero_points = 5;

leg = cell(1,length(sets_to_use));
instant_times = cell(1, length(sets_to_use));
instant_kr = cell(1, length(sets_to_use));
integ_kr = cell(1, length(sets_to_use));
max_instant_kr = 0;
max_integ_kr = 0;
idx = 1;
for set = sets_to_use
    times = [];
    instant_krs = [];
    integ_krs = [];
    set_struct = trace_struct_final([trace_struct_final.setID] == set);
    set_struct = set_struct([set_struct.MeanAP] > 1);
    for trace = set_struct
        tr_fluo = trace.fluo_interp;
        tr_times = trace.time_interp;
        tr_protein = trace.protein_interp;
        max_instant_kr = max([max_instant_kr trace.protein_interp]);
        nonzero_idxes = find(trace.fluo_interp > 0);
        if length(nonzero_idxes) > min_nonzero_points
            last_idx = length(nonzero_idxes);
            while last_idx >= min_burst_size && ...
                    nonzero_idxes(last_idx) - ...
                    nonzero_idxes(last_idx - min_burst_size + 1) >= min_burst_size
                last_idx = last_idx - 1;
            end
            if last_idx >= min_burst_size && last_idx ~= length(tr_fluo)
                times = [times tr_times(nonzero_idxes(last_idx))];
                instant_krs = [instant_krs tr_protein(nonzero_idxes(last_idx))];
                integ_kr = sum(tr_protein(1:nonzero_idxes(last_idx)));
                integ_krs = [integ_krs integ_kr];
                max_integ_kr = max(max_integ_kr, integ_kr);
            end
        end
    end
    instant_times{idx} = times;
    instant_kr{idx} = instant_krs;
    integ_kr{idx} = integ_krs;
    leg{idx} = ['set ' num2str(idx)];
    idx = idx + 1;
end

%scatter of instantaneous kr vs time
figure();
for i = 1:(idx - 1)
    scatter(instant_times{i}, instant_kr{i});
    hold on
end
title('Off times vs Instantaneous kr level');
xlabel('time');
ylabel('kruppel');
legend(leg);

% scatter of integrated kr vs time
figure();
for i = 1:(idx - 1)
    scatter(instant_times{i}, integ_kr{i});
    hold on
end
title('Off times vs Integrated kr level');
xlabel('time');
ylabel('kruppel');
legend(leg);

% histogram of off times
for i = 1:(idx-1)
    figure();
    histogram(instant_times{i});
    title(['Distribution of off times set ' num2str(sets_to_use(i))]);
end

% bins instantaneous and integrated kruppel levels


% script looks at instantaneous and integrated kruppel level when eve2
% turns off for the last time (eve2 is off once no more spots are detected
% for more than 4 time steps)

%imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

sets_to_use = [0]; %set 0 is all of them
all_sets = [1 4:6 8:12];
%all_sets = 1;
thresh = 0;
coeff = 2;
normalizer = 100;

leg = cell(1,length(sets_to_use));
instant_times = cell(1, length(sets_to_use));
instant_kr_cell = cell(1, length(sets_to_use));
integ_kr_cell = cell(1, length(sets_to_use));
idx = 1;
for set = sets_to_use
    times = [];
    instant_krs = [];
    integ_krs = [];
    if set == 0
        set_struct = trace_struct_final(ismember([trace_struct_final.setID], all_sets));
    else
        set_struct = trace_struct_final([trace_struct_final.setID] == set);
    end
    set_struct = set_struct([set_struct.MeanAP] > 1);
    for trace = set_struct
        tr_fluo = trace.fluo_interp;
        tr_times = trace.time_interp;
        tr_protein = trace.protein_interp / normalizer; % so it doesn't get too big
        tr_protein = tr_protein .^ coeff;
        above_thresh_idxes = find(tr_fluo > thresh);
        if isempty(above_thresh_idxes)
            continue
        end
        last_idx = above_thresh_idxes(end);
        if last_idx < length(tr_fluo)
            times = [times tr_times(last_idx)];
            instant_krs = [instant_krs tr_protein(last_idx)];
            integ_kr = sum(tr_protein(1:last_idx));
            integ_krs = [integ_krs integ_kr];
        end
    end
    instant_times{idx} = times;
    instant_kr_cell{idx} = instant_krs;
    integ_kr_cell{idx} = integ_krs;
    leg{idx} = ['set ' num2str(idx)];
    idx = idx + 1;
end

%scatter of instantaneous kr vs time
figure();
for i = 1:(idx - 1)
    scatter(instant_times{i}, instant_kr_cell{i});
    hold on
end
title('Off times vs Instantaneous kr level');
xlabel('time');
ylabel('kruppel');
legend(leg);

% scatter of integrated kr vs time
figure();
for i = 1:(idx - 1)
    scatter(instant_times{i}, integ_kr_cell{i});
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

for set = sets_to_use
    
    if set == 0
        set_struct = trace_struct_final(ismember([trace_struct_final.setID], all_sets));
    else
        set_struct = trace_struct_final([trace_struct_final.setID] == set);
    end
    set_struct = set_struct([set_struct.MeanAP] > 1);
    max_instant_kr = 0;
    max_integ_kr = 0;
    for trace = set_struct
        prot = (trace.protein_interp / normalizer) .^ coeff;
        max_instant_kr = max([max_instant_kr prot]);
        max_integ_kr = max(max_integ_kr, sum(prot));
    end
    max_instant_kr = max_instant_kr + 1;
    max_integ_kr = max_integ_kr + 1;

    % bins instantaneous and integrated kruppel levels
    num_bins = 50;
    instant_kr_jump = max_instant_kr / num_bins;
    instant_kr_divs = 0:instant_kr_jump:(max_instant_kr - .00001);
    instant_kr_counts = zeros(1, length(instant_kr_divs));
    instant_kr_bins = zeros(1, length(instant_kr_divs));
    integ_kr_jump = max_integ_kr / num_bins;
    integ_kr_divs = 0:integ_kr_jump:(max_integ_kr-.00001);
    integ_kr_counts = zeros(1, length(integ_kr_divs));
    integ_kr_bins = zeros(1, length(integ_kr_divs));
    
    for trace = set_struct
        tr_protein = trace.protein_interp / normalizer;
        tr_protein = tr_protein .^ coeff;
        tr_fluo = trace.fluo_interp;
        above_thresh = find(tr_fluo > thresh);
        if isempty(above_thresh)
            continue
        end
        last_idx = above_thresh(end);
        if last_idx == length(trace.protein_interp)
            continue
        end
        integ_so_far = 0;
        for i = 1:last_idx
            instant_kr_bin = floor(tr_protein(i) / instant_kr_jump) + 1;
            instant_kr_counts(instant_kr_bin) = instant_kr_counts(instant_kr_bin) + 1;
            integ_so_far = integ_so_far + tr_protein(i);
            integ_kr_bin = floor(integ_so_far / integ_kr_jump) + 1;
            integ_kr_counts(integ_kr_bin) = integ_kr_counts(integ_kr_bin) + 1;
        end
        if last_idx < length(tr_protein)
            instant_kr_bins(instant_kr_bin) = instant_kr_bins(instant_kr_bin) + 1;
            integ_kr_bins(integ_kr_bin) = integ_kr_bins(integ_kr_bin) + 1;
        end
    end
    integ_kr_cum = cumsum(integ_kr_bins);
    integ_kr_cum = integ_kr_cum / sum(integ_kr_bins);
    integ_kr_prob = integ_kr_cum(1:end) - [0 integ_kr_cum(1:end-1)];

    instant_kr_counts(instant_kr_counts == 0) = 1;
    integ_kr_counts(integ_kr_counts == 0) = 1;
    instant_kr_bins = instant_kr_bins ./ instant_kr_counts;
    integ_kr_bins = integ_kr_bins ./ integ_kr_counts;
    
    figure();
    plot(instant_kr_divs, instant_kr_bins);
    title(['Instant Kruppel Conc turn off, set ' num2str(set)]);
    figure();
    plot(integ_kr_divs, integ_kr_bins);
    title(['Integrated Kruppel Conc turn off, set ' num2str(set)]);
    
    figure();
    plot(instant_kr_divs, instant_kr_counts);
    title(['Instant Kruppel Conc counts, set ' num2str(set)]);
    figure();
    plot(integ_kr_divs, integ_kr_counts);
    title(['Integrated Kruppel Conc counts, set ' num2str(set)]);
    
    figure();
    plot(integ_kr_divs, integ_kr_cum);
    title(['Integrated Kruppel Conc turn off cum, set ' num2str(set)]);
    figure();
    plot(integ_kr_divs, integ_kr_prob);
    title(['Integrated Kruppel Conc turn off prob, set ' num2str(set)]);
    
    
end

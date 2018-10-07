% compares standard deviation by ap at a specific time

% import fixed data
raw = csvread('../dat/Kruppel_fixed/Raw_Profiles/Data_Kr_raw.csv',3,0);
proc = csvread('../dat/Kruppel_fixed/Processed_Profiles/Data_Kr_proc.csv',4,0);

% import live data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';

sets_to_use = [2 3 4 5 9 12];
ap_bins_live = -3:0.5:10;
ap_bins_fixed = 365:5:495;
live_width = 0.5;
fixed_width = 5;
min_time = 27 * 60;
max_time = 30 * 60;

krup_live = cell(1, length(ap_bins_live));
krup_fixed = cell(1, length(ap_bins_fixed));

% collects all relevant kruppel points for each ap bin for live data
set_struct = trace_struct_final(ismember([trace_struct_final.setID], sets_to_use));
set_struct = set_struct([set_struct.MeanAP] > -3);
set_struct = set_struct([set_struct.MeanAP] < 10.5);
for trace = set_struct
    idxes = find(trace.time_interp >= min_time & trace.time_interp <= max_time);
    if isempty(idxes)
        continue;
    end
    ap_idx = floor((trace.MeanAP + 3) / live_width) + 1;
    krup_live{ap_idx} = [krup_live{ap_idx} trace.protein_interp(idxes)];
end

% collects all relevant kruppel points for fixed data
for i = 1:size(raw,1)
    time = raw(i,4);
    sub = min(raw(i,5:end));
    if isnan(time) || time >= max_time / 60 || time < min_time / 60
        continue
    end
    for j = 1:length(ap_bins_fixed)
        ap_start = ap_bins_fixed(j);
        val = sum(raw(i, ap_start:(ap_start + fixed_width - 1)));
        krup_fixed{j} = [krup_fixed{j} (val - sub * fixed_width)];
    end 
end

% consolidates data
krup_live_stds = zeros(1, length(krup_live));
krup_fixed_stds = zeros(1, length(krup_fixed));
for i = 1:length(krup_live)
    krup_live_stds(i) = std(krup_live{i}) ./ mean(krup_live{i});
    krup_fixed_stds(i) = std(krup_fixed{i}) ./ mean(krup_fixed{i});
end

% plots data
figure();
plot(ap_bins_live, krup_live_stds);
hold on
plot(ap_bins_live, krup_fixed_stds);
ylabel('kruppel standard deviation');
xlabel('ap bin');
legend('live data', 'fixed data');
title('Standard deviation of heterozygous kruppel variation 27-30 minutes');

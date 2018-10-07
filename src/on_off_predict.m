% performs a multivariate linear regression to predict average on/off
% durations as a function of ap, time, and kruppel

% imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

min_ap = -0.5;
trace_struct_final = trace_struct_final([trace_struct_final.MeanAP] > min_ap);
sets_to_use = [1 4:6 8:12];

half_elong = 3;
thresh = 50000;
num_bins = 20;

% sets up bins for krupel, time, and ap
all_krup = [trace_struct_final.protein_interp];
all_krup = all_krup(all_krup > 0);
min_krup = min(all_krup);
max_krup = max(all_krup);
krup_step = (max_krup - min_krup) / num_bins;
krup_bins = min_krup:krup_step:(max_krup - 1);

all_times = [trace_struct_final.time_interp];
min_time = 0;
max_time = max(all_times);
time_step = (max_time - min_time) / num_bins;
time_bins = min_time:time_step:(max_time - 1);

all_ap = [trace_struct_final.MeanAP];
min_ap = min(all_ap);
max_ap = max(all_ap);
ap_step = (max_ap - min_ap) / num_bins;
ap_bins = min_ap:ap_step:(max_ap - 1);

all_on = cell([length(krup_bins), length(time_bins), length(ap_bins)]);
all_off = cell([length(krup_bins), length(time_bins), length(ap_bins)]);

% bins values from each trace for later averaging
set_struct = trace_struct_final(ismember([trace_struct_final.setID], sets_to_use));
for trace = set_struct
    ap_val = trace.MeanAP;
    ap_idx = floor((ap_val - min_ap) / ap_step) + 1;
    ap_idx = max(min(length(ap_bins), ap_idx),1);
    below_thresh = find(trace.fluo_interp < thresh);
    for i = 1:(length(below_thresh) - 1)
        if below_thresh(i+1) - below_thresh(i) > 1
            krup_val = mean(trace.protein_interp(...
                below_thresh(i):below_thresh(i+1)));
            krup_idx = floor((krup_val - min_krup) / krup_step) + 1;
            krup_idx = max(min(length(krup_bins), krup_idx),1);
            time_val = trace.time_interp(below_thresh(i+1));
            time_idx = floor((time_val - min_time) / time_step) + 1;
            time_idx = max(min(length(time_bins), time_idx),1);
            
            all_on{krup_idx, time_idx, ap_idx} = [all_on{krup_idx,time_idx,ap_idx} ...
                below_thresh(i+1) - below_thresh(i) - 1];
        end
    end
    above_thresh = find(trace.fluo_interp > thresh);
    for i = 1:(length(above_thresh) - 1)
        if above_thresh(i+1) - above_thresh(i) > 1
            krup_val = mean(trace.protein_interp(...
                above_thresh(i):above_thresh(i+1)));
            krup_idx = floor((krup_val - min_krup) / krup_step) + 1;
            krup_idx = max(min(length(krup_bins), krup_idx),1);
            time_val = trace.time_interp(above_thresh(i+1));
            time_idx = floor((time_val - min_time) / time_step) + 1;
            time_idx = max(min(length(time_bins), time_idx),1);
            all_off{krup_idx, time_idx, ap_idx} = [all_off{krup_idx,time_idx,ap_idx} ...
                above_thresh(i+1) - above_thresh(i) - 1];
        end
    end
end

% put data in usable form for regression
X = [];
Y_on = [];
Y_off = [];
min_points = 5;
idx = 1;
for i = 1:length(krup_bins)
    for j = 1:length(time_bins)
        for k = 1:length(ap_bins)
            on_points = all_on{i,j,k};
            off_points = all_off{i,j,k};
            if length(on_points) >= min_points
                X = [X; krup_bins(i), 0, ap_bins(k), 1];
                Y_on(idx) = mean(on_points);
                Y_off(idx) = mean(off_points);
                idx = idx + 1;
            end
        end
    end
end

[beta_on,Sigma_on,E_on,CovB_on,logL_on] = mvregress(X, Y_on');
[beta_off,Sigma_off,E_off,CovB_off,logL_off] = mvregress(X, Y_off');

% looks at the variance in raw and processed fixed kruppel data

%import data
raw = csvread('../dat/Kruppel_fixed/Raw_Profiles/Data_Kr_raw.csv',3,0);
proc = csvread('../dat/Kruppel_fixed/Processed_Profiles/Data_Kr_proc.csv',4,0);

% plots mean/std of raw data at 3 minute intervals as function of
% time for various ap positions
ap_rad = 0;
ap_bins = 365:30:485;
time_bins = 0:3:48;
leg = {};
idx = 1;
ap_means = cell(1, length(ap_bins));
ap_stds = cell(1, length(ap_bins));
ap_means_sub = cell(1, length(ap_bins));
ap_stds_sub = cell(1, length(ap_bins));
ap_maxes_sub = cell(1, length(ap_bins));
ap_mins_sub = cell(1, length(ap_bins));
ap_means_scaled = cell(1, length(ap_bins));
ap_stds_scaled = cell(1, length(ap_bins));
scatter_times = cell(1, length(ap_bins));
scatter_vals_sub = cell(1, length(ap_bins));
ap_means1 = cell(1, length(ap_bins));
ap_stds1 = cell(1, length(ap_bins));
ap_means2 = cell(1, length(ap_bins));
ap_stds2 = cell(1, length(ap_bins));
figure();
for ap = ap_bins
    fluos = cell(1, length(time_bins));
    fluos_sub = cell(1, length(time_bins));
    fluos_scaled = cell(1, length(time_bins));
    fluos1 = cell(1, length(time_bins));
    fluos2 = cell(1, length(time_bins));
    counts = zeros(1, length(time_bins));
    for i = 1:size(raw,1)
        time = raw(i,4);
        sub = min(raw(i,5:end));
        scale = max(raw(i,5:end));
        if isnan(time) || time >= (time_bins(end) + 3)
            continue
        end
        time_idx = floor(time / 3) + 1;
        for j = ap - ap_rad:ap + ap_rad
            fluos{time_idx} = [fluos{time_idx} raw(i, 5 + j)];
            fluos_sub{time_idx} = [fluos_sub{time_idx} raw(i, 5 + j) - sub];
            fluos_scaled{time_idx} = [fluos_scaled{time_idx} ...
                (raw(i, 5 + j) - sub) / scale];
            counts(time_idx) = counts(time_idx) + 1;
            scatter_times{idx} = [scatter_times{idx} time];
            scatter_vals_sub{idx} = [scatter_vals_sub{idx} raw(i, 5 + j) - sub];
            if raw(i,2) == 1
                fluos1{time_idx} = [fluos1{time_idx} raw(i, 5+j) - sub];
            else
                fluos2{time_idx} = [fluos2{time_idx} raw(i,5+j) - sub];
            end
        end 
    end
    for j = 1:length(fluos)
        if ~isempty(fluos{j})
            ap_means{idx}(j) = mean(fluos{j});
            ap_stds{idx}(j) = std(fluos{j});
            ap_means_sub{idx}(j) = mean(fluos_sub{j});
            ap_stds_sub{idx}(j) = std(fluos_sub{j});
            ap_maxes_sub{idx}(j) = max(fluos_sub{j});
            ap_mins_sub{idx}(j) = min(fluos_sub{j});
            ap_means_scaled{idx}(j) = mean(fluos_scaled{j});
            ap_stds_scaled{idx}(j) = std(fluos_scaled{j});
            ap_means1{idx}(j) = mean(fluos1{j});
            ap_stds1{idx}(j) = std(fluos1{j});
            ap_means2{idx}(j) = mean(fluos2{j});
            ap_stds2{idx}(j) = std(fluos2{j});
        end
    end
%     display(['ap ' num2str(ap)]);
%     rel_std = sort(ap_stds_sub{idx} ./ ap_means_sub{idx});
%     display(rel_std);
%     rel_range = sort((ap_maxes_sub{idx} - ap_mins_sub{idx}) ./ ap_means_sub{idx});
%     display(rel_range);
    leg{idx} = num2str(ap / 10);
    idx = idx + 1;
    plot(time_bins, counts);
    hold on
end
title('raw counts per time');
legend(leg);

% plots non adjusted
figure();
for idx = 1:length(ap_bins)
    errorbar(time_bins, ap_means{idx}, ap_stds{idx});
    hold on
end
title('Kruppel content non adjusted');
legend(leg);
xlabel('time (min)');
ylabel('a.u.');

% plots background subtracted
figure();
for idx = 1:length(ap_bins)
    errorbar(time_bins, ap_means_sub{idx}, ap_stds_sub{idx});
    hold on
end
title('Kruppel content background subtracted');
legend(leg);
xlabel('time (min)');
ylabel('a.u.');

% plots scaled
figure();
for idx = 1:length(ap_bins)
    errorbar(time_bins, ap_means_scaled{idx}, ap_stds_scaled{idx});
    hold on
end
title('Kruppel content scaled');
legend(leg);
xlabel('time (min)');
ylabel('a.u.');

% plots scatter of background subtracted values
for idx = 1:length(ap_bins)
    figure();
    scatter(scatter_times{idx}, scatter_vals_sub{idx});
    title(['scatter of kruppel vals vs time ap: ' leg{idx}]);
    xlabel('time (min)');
    ylabel('a.u.');
end

% plot 1 and 2 separately for background subtracted kruppel

figure();
for idx = 1:length(ap_bins)
    errorbar(time_bins, ap_means1{idx}, ap_stds1{idx});
    hold on
end
title('Kruppel content background subtracted set 1');
legend(leg);
xlabel('time (min)');
ylabel('a.u.');

figure();
for idx = 1:length(ap_bins)
    errorbar(time_bins, ap_means2{idx}, ap_stds2{idx});
    hold on
end
title('Kruppel content background subtracted set 2');
legend(leg);
xlabel('time (min)');
ylabel('a.u.');



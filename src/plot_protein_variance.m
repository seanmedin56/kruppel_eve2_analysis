% looks at nearest ap neighbor variance for each set (as function of mean
% fluorescence between the two neighbors and also plots the distribution of 
% jumps in kruppel level from time step to time step

% loads data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);
sets_to_use = 1:3;
figure();
for set = sets_to_use
    
    trace_struct_set = trace_struct_final([trace_struct_final.setID] == set);
    
    % find maximum and minimum kruppel (for binning)
    num_bins = 30;
    max_krup = 0;
    min_krup = 1000000;
    for trace = trace_struct_final
        max_krup = max([trace.protein_interp max_krup]);
        min_krup = min([trace.protein_interp, min_krup]);
    end

    bin_size = (max_krup - min_krup) / num_bins;
    krup_bins = min_krup:bin_size:(max_krup - bin_size);

    % finds nearest ap neighbor for each trace and finds differences in kruppel
    % level at each time between them
    krups_cell = cell(1, length(krup_bins));

    for trace = trace_struct_set
        ap = trace.MeanAP;
        [M,I] = sort(abs([trace_struct_set.MeanAP] - ap));
        closest_trace = trace_struct_set(I(2));
        [times_interp,ia,ib] = intersect(trace.time_interp, closest_trace.time_interp);
        new_dists = abs(trace.protein_interp(ia) - closest_trace.protein_interp(ib));
        new_avgs = (trace.protein_interp(ia) + closest_trace.protein_interp(ib)) * .5;
        for i = 1:length(new_dists)
            dist = new_dists(i);
            avg = new_avgs(i);
            krups_cell{floor((avg - min_krup) / bin_size + 1)} = ... 
                [krups_cell{floor((avg - min_krup) / bin_size + 1)} dist];
        end  
    end

    krup_means = zeros(1,length(krups_cell));
    krup_stds = zeros(1,length(krups_cell));
    for i = 1:length(krups_cell)
        krup_means(i) = mean(krups_cell{i});
        krup_stds(i) = std(krups_cell{i});
    end
    errorbar(krup_bins, krup_means, krup_stds);
    hold on
    
end

title('Kruppel nearest neighbor distance vs kruppel content');
xlabel('mean fluorescence');
ylabel('mean fluorescence distance');

% plots histogram of distance between adjacent time points for each set
for set = sets_to_use
    adj_points = [];
    trace_struct_set = trace_struct_final([trace_struct_final.setID] == set);
    for trace = trace_struct_set
        adj_points = [adj_points diff(trace.protein_interp)];
    end
    figure();
    histogram(adj_points);
    title(['Adjacent Times Fluo Distribution Set ' num2str(set)]);
    display(mean(adj_points));
    display(std(adj_points));
end


% takes middle parts of each trace and clusters them

%imports data

time_min = 600;
time_max = 2100;
sets_to_use = [1 4:6 8:12];
num_clusters = 10;
traces_mat = [];
idx = 1;
% formats traces for clustering
for set = sets_to_use
    set_struct = trace_struct_final([trace_struct_final.setID] == set);
    for trace = set_struct
        if ismember(time_min, trace.time_interp) && ...
                ismember(time_max, trace.time_interp)
            start_idx = time_min == trace.time_interp;
            end_idx = time_max == trace.time_interp;
            traces_mat(idx) = trace.protein_interp(start_idx:end_idx); 
            idx = idx + 1;
        end
    end
end

clusters = kmeans(traces_mat, num_clusters);

cluster_struct = struct;
for i = 1:num_clusters
    
end

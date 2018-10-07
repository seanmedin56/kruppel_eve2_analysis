% tries to find how much transcription there is at any given time point by
% only assuming a constant elongation time

%imports data
data_path = '../dat/Kruppel_eve2_pass1/inference_traces_Kruppel_eve2_pass1_dT20.mat';
load(data_path);

sets_to_use = [8];
elong_time = 8; %number of time points before fall off
rise_time = 2;
buffer = 8;

set_struct = trace_struct_final(ismember([trace_struct_final.setID], sets_to_use));
all_points = [];
for trace = set_struct
    num_points = length(trace.time_interp);
    A = zeros(num_points,num_points - buffer);
    for i = 1:(num_points - buffer)
        for j = 1:elong_time
            if i + j - 1 > num_points
                break
            end
            if j <= rise_time
                A(i + j - 1, i) = j /( 2 * rise_time);
            else
                A(i + j - 1, i) = 1;
            end
        end
    end
    fluo_tr = trace.fluo_interp';
    x = lsqnonneg(A, fluo_tr);   
    nonzero_idxes = find(x > 0);
    all_points = [all_points x(nonzero_idxes)'];
end
figure();
histogram(all_points);
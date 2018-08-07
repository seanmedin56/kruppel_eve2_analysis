function [lsq_fun] = predict_eve2(gt, kr, bc, hb, gt_ks, kr_ks, bc_ks, ...
    hb_ks, time_res, elong_window, fluo_per_time)
%predict_eve2: calculates expected amount of fluorescence from a nuceli 
% from eve2 transcription for each time step
%   Assumes constant protein content over single time step
%   Agruments: gt,kr,bc,hb: protein concentration at each time step
%              gt_ks,kr_ks,bc_ks,hb_ks: reaction coefficients for each
%              protein

    lsq_fun = zeros(1,length(gt));
    cur_eve2 = zeros(1, elong_window);
    cur_prob = zeros(1,16);
    cur_prob(1) = 1; % starts out with no proteins bound to promoter
    for i = 1:length(gt)
        cur_eve2 = [cur_eve2(2:end) cur_prob(5) * fluo_per_time];
        cur_prob = calc_probs(gt(i), kr(i), bc(i), hb(i), gt_ks, kr_ks, ...
            bc_ks, hb_ks, time_res, cur_prob);
        lsq_fun(i) = sum(cur_eve2);
    end

end


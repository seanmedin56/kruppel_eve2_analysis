function [probs] = calc_probs(gt, kr, bc, hb, gt_ks, kr_ks, bc_ks, ...
    hb_ks, time_res, prior_probs)
%calc_probs: calculates probabilities fo being in a certain state for eve2
%regulatory system
%   Assumes constant protein content over time step
%   prob = [p(), p(Hb), p(Bc), p(HbBc), p(Kr), p(HbKr), p(BcKr), ...]
%   Agruments: gt,kr,bc,hb: protein concentration at the time step
%              gt_ks,kr_ks,bc_ks,hb_ks: reaction coefficients for each
%              protein
          
    R = zeros(16,16);
    hb_idxes = 2:2:16;
    no_hb_idxes = 1:2:16;
    bc_idxes = [3:4:16 4:4:16];
    no_bc_idxes = [1:4:16 2:4:16];
    kr_idxes = [5:8 13:16];
    no_kr_idxes = [1:4 9:12];
    gt_idxes = 9:16;
    no_gt_idxes = 1:8;
    R(no_hb_idxes, hb_idxes) = R(no_hb_idxes, hb_idxes) + hb_ks(2);
    R(hb_idxes, no_hb_idxes) = R(hb_idxes, no_hb_idxes) + hb_ks(1) * hb;
    R(no_bc_idxes, bc_idxes) = R(no_bc_idxes, bc_idxes) + bc_ks(2);
    R(bc_idxes, no_bc_idxes) = R(bc_idxes, no_bc_idxes) + bc_ks(1) * bc;
    R(no_kr_idxes, kr_idxes) = R(no_kr_idxes, kr_idxes) + kr_ks(2);
    R(kr_idxes, no_kr_idxes) = R(kr_idxes, no_kr_idxes) + kr_ks(1) * kr;
    R(no_gt_idxes, gt_idxes) = R(no_gt_idxes, gt_idxes) + gt_ks(2);
    R(gt_idxes, no_gt_idxes) = R(gt_idxes, no_gt_idxes) + gt_ks(1) * gt;
    probs = expm(R * time_res) * prior_probs;
end


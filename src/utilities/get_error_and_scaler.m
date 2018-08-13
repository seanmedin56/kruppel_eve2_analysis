function [error, scaler] = get_erorr_and_scaler(fit_mat, other_mats, ...
                  other_time_delays, other_ap_delays, ap_delay, time_delay)
%GET_ERORR_AND_SCALER Summary of this function goes here
%   Detailed explanation goes here

    % finds best scaler
    upper = 0;
    lower = 0;
    for i = 1:size(fit_mat,1)
        for j = 1:size(fit_mat,2)
            lower_contrib = 0;
            for k = 1:length(other_ap_delays)
                ap_other = ap_delay - other_ap_delays(k) + i;
                time_other = time_delay - other_time_delays(k) + j;
                if ap_other > 0 && time_other > 0 && ...
                        ap_other < size(other_mats{k},1) && ... 
                        time_other < size(other_mats{k},2)
                    lower_contrib = lower_contrib + 1;
                    upper = upper + other_mats{k}(ap_other,time_other) * ...
                        fit_mat(i,j);
                end
            end
            lower = lower + lower_contrib * fit_mat(i,j)^2;
            if isnan(lower)
                stop = true;
            end
        end
    end
    scaler = upper / lower;
    
    % calculates error
    counts = 0;
    error = 0;
    for i = 1:size(fit_mat,1)
        for j = 1:size(fit_mat,2)
            for k = 1:length(other_ap_delays)
                ap_other = ap_delay - other_ap_delays(k) + i;
                time_other = time_delay - other_time_delays(k) + j;
                if ap_other > 0 && time_other > 0 && ...
                        ap_other < size(other_mats{k},1) && ... 
                        time_other < size(other_mats{k},2)
                    error = error + (other_mats{k}(ap_other,time_other) ...
                        - fit_mat(i,j))^2;
                    counts = counts + 1;
                    if isnan(error)
                        stop = true;
                    end
                end
            end
        end
    end
    error = error / counts / counts;
    if isnan(error)
        stop = true;
    end
end


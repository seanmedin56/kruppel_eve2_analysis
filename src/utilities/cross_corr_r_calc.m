function cross_corr_r = cross_corr_r_calc(trace1, trace2, max_delay)

% returns the average raw moment crosscorrelation of the traces

% ----------calculates individual correaltions (with weights)--------------

    corrs = cell([1 length(trace1)]);
    for i = 1:length(trace1)
        limit = min([max_delay, length(trace1{i}), length(trace2{i})]);
        len = min(length(trace1{i}), length(trace2{i}));
        corr = zeros([1 2*limit - 1]);
        for j = 1:limit
            corr(j) = trace1{i}(limit - j + 1:len) * trace2{i}(1:len-limit+j)';
        end
        for j = 2:limit
            corr(j + limit - 1) = trace1{i}(1:len - j + 1) * trace2{i}(j:len)';
        end
        corr = corr ./ [(len - limit + 1):len (len - 1:-1:(len - limit + 1))];
        corrs{i} = corr / max(abs(corr));
    end
    
% -----------------combines correaltions together------------------------    

    cross_corr_r = zeros([1 2 * limit - 1]);
    for i = 1:length(corrs)
        cross_corr_r(1:length(corrs{i})) = cross_corr_r(1:length(corrs{i})) ...
            + corrs{i};
    end
    cross_corr_r = cross_corr_r / max(abs(cross_corr_r));
end


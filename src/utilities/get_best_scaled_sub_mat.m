function [best_mat,best_scalar,row,col] = get_best_scaled_sub_mat(A_small, B_large)
%GET_BEST_SCALED_SUB_MAT Summary of this function goes here
%   Detailed explanation goes here
    
    max_row_start = size(B_large,1) - size(A_small,1) + 1;
    max_column_start = size(B_large,2) - size(A_small,2) + 1;
    best_mat = zeros(size(A_small));
    best_scalar = 0;
    row = 1;
    col = 1;
    best_error = norm(A_small - best_mat, 'fro');
    for i = 1:max_row_start
        for j = 1:max_column_start
            test_mat = B_large(i:i+size(A_small,1)-1, j:j+size(A_small,2)-1);
            scalar = sum(sum(test_mat .* A_small)) / sum(sum(test_mat .* test_mat));
            error = norm(test_mat * scalar - A_small, 'fro');
            if best_error > error
                best_error = error;
                best_mat = test_mat * scalar;
                best_scalar = scalar;
                row = i;
                col = j;
            end
        end
    end

end


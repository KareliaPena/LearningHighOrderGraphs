function [A]= symmetrize_tensor(A)
    %%%take a p order tensor (n x n x n_3 x ...... x n_p) as input, modify it to be its symmetric version with
    %%%dimension (n x n x (2n_3 + 1) x ...... x (2n_p +1)). The symmetric() operation should be applied to each 
    %%%order of the tensor.
    
    p = length(size(A));
    for i=3:p
        A = sym_order(A, i); %apply to each dimension from 3 to np
    end
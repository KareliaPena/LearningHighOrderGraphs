function At = t_tran(A)
    %%%take a p order tensor (n x n x n_3 x ...... x n_p) as input, compute the transpose.
    %%The tran() operation should be applied to each order of the tensor.%%%

    if length(size(A))==3
       At = tran3(A);
       return
    end
    sz = size(A);
    sz2= sz;
    sz2(1)=sz(2);
    sz2(2)=sz(1);
    At=zeros(sz2);
    for i= 1:sz(end)
        otherdims = repmat({':'},1,ndims(A)-1);
        At(otherdims{:},i) = t_tran(A(otherdims{:},i)); %apply to each dimension from 3 to np
    end
end
function A_order=sym_order(A, i)
    %%%"Apply the symmetric operation to the ith dimension.%%%
    Ashape=size(A);
    first_0 = zeros([Ashape(1:i-1),1,Ashape(i+1:end)]);
    reverse_A = flip(double(A), i);
    A_order = 1/2*sptensor(cat(i,first_0,double(A), double(reverse_A)));
end
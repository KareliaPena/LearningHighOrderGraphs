
function [n_R] = normR(R)

%R = full(R);

M_R = R*R';

a = M_R(1);
b = M_R(2);
N = size(R,1);

n_R = (a+b*(N-1))^(1/2);

end

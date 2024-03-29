function Xt = tran3(X)

% conjugate transpose of a 3 way tensor 
% X  - n1*n2*n3 tensor
% Xt - n2*n1*n3  tensor
%
% version 1.0 - 18/06/2016
%
% Written by Canyi Lu (canyilu@gmail.com)
% 
%
% References: 
% Canyi Lu, Tensor-Tensor Product Toolbox. Carnegie Mellon University. 
% June, 2018. https://github.com/canyilu/tproduct.
%
% Canyi Lu, Jiashi Feng, Yudong Chen, Wei Liu, Zhouchen Lin and Shuicheng
% Yan, Tensor Robust Principal Component Analysis with A New Tensor Nuclear
% Norm, arXiv preprint arXiv:1804.03728, 2018
%

sz = size(X);
Xt = zeros(sz(2),sz(1),sz(3));
Xt(:,:,1) = X(:,:,1)';
for i = 2 : sz(3)
    Xt(:,:,i) = X(:,:,sz(3)-i+2)';
end
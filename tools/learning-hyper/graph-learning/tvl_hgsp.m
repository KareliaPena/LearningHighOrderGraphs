function [L_res,L,costFinal_MAP,costFinal_MSE] = tvl_hgsp(Y,alpha,beta)

%% Laplacian constraints
M = length(size(Y));
N = size(Y,1);

try
    load(['L_N=' num2str(N)]);

catch
    laplacian_constraints_teig_spatial_non_uniform(N,M)
    load(['L_N=' num2str(N)]);
end

% mat_obj -> Operator that accounts for symmetric: from vector to whole tensor
% mat_obj2 -> This math object accounts for the symmetry in the third dimenssion
% mat_obj3 -> This math object maps the data cube organized from frontal slices to scalar tubes
% mat_obj4  -> This math object maps the data cube organized from scalar tubes to frontal slices
% mat_obj_fft -> maps the data cube organized in scalar tubes from the spatial domain to the fft domain
Y=Y/max(Y(:));
allp=[];
for i=1:2*N+1
    p = vec(squeeze(Y(:,:,i))*squeeze(Y(:,:,i))')';
    allp=[allp,p];
end
allp=reshape(allp,[N^2,2*N+1]);
% allp=allp/max(allp(:));
%% optimization
 cvx_begin
 
% cvx_solver mosek
 
variable L(NDE,1)  %L is in the spatial domain

L_tubes_tensor=mat_obj3*mat_obj2*mat_obj*L;
L_tubes_tensor=reshape(L_tubes_tensor,[size(L_tubes_tensor,1)/N^2,N^2]);
L_tubes_tensor_fft=mat_obj_fft*L_tubes_tensor;

L_fft=mat_obj4*L_tubes_tensor_fft(:);
%L_fft=L_fft(1:N^2*(N+1)); %keeps only different slices that is first N+1
L_fft=reshape(L_fft,[N^2,2*N+1]);
FinalTube_fft=diag(L_fft'*allp);
ainv = conj(mat_obj_fft)/(2*N+1);
FinalTube=ainv*FinalTube_fft;

minimize alpha*sum(FinalTube) + beta*sum_square_abs(mat_obj2*mat_obj*L)

subject to
    A1*L == b1
    A2*L <= b2
cvx_end

L_tubes_tensor=mat_obj3*mat_obj2*mat_obj*L;
L_tubes_tensor=reshape(L_tubes_tensor,[size(L_tubes_tensor,1)/N^2,N^2]);
L_tubes_tensor_fft=mat_obj_fft*L_tubes_tensor;

L_fft=mat_obj4*L_tubes_tensor_fft(:);
%L_fft=L_fft(1:N^2*(N+1)); %keeps only different slices that is first N+1
L_fft=reshape(L_fft,[N^2,2*N+1]);
FinalTube_fft=diag(L_fft'*allp);
ainv = conj(mat_obj_fft)/(2*N+1);
FinalTube=ainv*FinalTube_fft;

costFinal_MAP=sum(FinalTube);
costFinal_MSE=sum_square_abs(L);

%% convert from vector form to tensor form
L_res = reshape(mat_obj2*mat_obj*L,[N,N,2*N+1]);

function  [V,S,Shift_Operator_Abs_Norm,Shift_Operator_Norm] = teig3(A,ordereig)
%  t_eig_Asc   Find the t-eigendecomposition of matrix Yr
%
%  [Shift_Operator_Abs_Norm,Shift_Operator_Norm,V,S] = t_eig_Asc(A) returns a f-diagonal tensor S
%  of eigen-tuples and a tensor V of eigen-matrices.
%  Shift_Operator_Abs_Norm has positive eigentuples.
%  Shift_Operator_Norm has normalized eigentuples
%  dir determines the order of the eigen-tuples

SzT = size(A);
A = fft(double(A),[],3);
powerIm=sum(sum(sum(imag(A).^2)));
if powerIm<eps
    A=real(A);
end


V = zeros(SzT);
S = zeros(SzT);
Yt = zeros(SzT);
Yt2 = zeros(SzT);

for i = 1 : SzT(3)
    [V(:,:,i),S(:,:,i)]= eigs(A(:,:,i),SzT(1),ordereig);
    maxEig=max(max(abs(S(:,:,i))));
    Yt(:,:,i) = V(:,:,i)*(abs(S(:,:,i))/maxEig)*V(:,:,i)';
    Yt2(:,:,i) = V(:,:,i)*((S(:,:,i))/maxEig)*V(:,:,i)';     
end
Shift_Operator_Abs_Norm=ifft(Yt,[],3,'symmetric');
Shift_Operator_Norm=ifft(Yt2,[],3,'symmetric');

V = ifft(V,[],3,'symmetric');
S = ifft(S,[],3,'symmetric');
end
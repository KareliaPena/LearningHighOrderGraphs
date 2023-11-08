function [Zs,Zs_fft] = hgsp_pointwise_aggregate(A,B)
    % Take a tensors of order p, output their t eigendecomposition. 
    % Suppose A is in dimension n1 x n2 x ... x np.
    % ordereig = Largest Magnitude 'LM' or Smallest Magnitude 'SM' 
    A=double(A);
    tol=eps;
    shape_ten=size(A);
    p = length(shape_ten); % tensor order
    num_slices = prod(shape_ten(3:end)); % n3 x n4 x ... x np
    flatten_shape_A = [shape_ten(1:2), num_slices];
    flatten_shape_Z = [1,1, num_slices];
    shape_z = [1,1, shape_ten(3:end)];
    A_fft=double(A);
    B_fft=double(B);
    %%% Conduct fft along every dimension after 2 recursively
    for i=3:p %skip the first two dimension
        A_fft = fft(A_fft,[],i);
        B_fft = fft(B_fft,[],i);
    end

     %%% Unfold all dimensions after 2
    flatten_A = reshape(A_fft,flatten_shape_A);
    flatten_B = reshape(B_fft,flatten_shape_A);
    
    powerIm=sum(imag(flatten_A(:)).^2);
    if powerIm<tol
        flatten_A=real(flatten_A);
    end
    
    flatten_Z_fft = zeros(flatten_shape_Z);
    

    %%% Regular eigen value decomposition for each slices
    for i=1:num_slices
        rng(0)
        flatten_Z_fft(:,:,i) = sum(sum(flatten_A(:,:,i).*flatten_B(:,:,i)));
    end
    %%% Unfold all dimensions after 2
    Zs_fft = reshape(flatten_Z_fft,shape_z);
    Zs=Zs_fft;
    for i=p:-1:3 %in reversed(range(2, p)):
        Zs = ifft(Zs,[], i);
    end
    powerIm=sum(imag(Zs(:)).^2);
    if powerIm<tol
        Zs=real(Zs);
    end
    
end
    
    
    
    
    
  

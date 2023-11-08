function [Zs,Zs_fft] = hgsp_distanz(A)
    % Take a tensors of order p, output their t eigendecomposition. 
    % Suppose A is in dimension n1 x n2 x ... x np.
    % ordereig = Largest Magnitude 'LM' or Smallest Magnitude 'SM' 
    A=double(A);
    tol=eps;
    shape_ten=size(A);
    p = length(shape_ten); % tensor order
    num_slices = prod(shape_ten(3:end)); % n3 x n4 x ... x np
    flatten_shape_D = [shape_ten(1:2), num_slices];
    
    shape_z=shape_ten;
    shape_z(2)=shape_z(1);
    flatten_shape_Z = [shape_z(1:2), num_slices];
    D=double(A);
    %%% Conduct fft along every dimension after 2 recursively
    for i=3:p %skip the first two dimension
        D = fft(D,[],i);
    end

     %%% Unfold all dimensions after 2
    flatten_D = reshape(D,flatten_shape_D);
  
    powerIm=sum(imag(flatten_D(:)).^2);
    if powerIm<tol
        flatten_D=real(flatten_D);
    end
    
    flatten_Z_fft = zeros(flatten_shape_Z);
    

    %%% Regular eigen value decomposition for each slices
    for i=1:num_slices
        rng(0)
        flatten_Z_fft(:,:,i) = gsp_distanz(flatten_D(:,:,i)').^2;
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
    
    
    
    
    
  

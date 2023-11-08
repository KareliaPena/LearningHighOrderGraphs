function [flatten_D,D] = t_fft(A)
    % Computes the fft on diffeent directions
    tol=eps;
    shape_ten=size(A);
    p = length(shape_ten); % tensor order
    num_slices = prod(shape_ten(3:end)); % n3 x n4 x ... x np
    flatten_shape_D = [shape_ten(1:2), num_slices];
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
        D=real(D);
    end
  

    
    
    
    
    
  

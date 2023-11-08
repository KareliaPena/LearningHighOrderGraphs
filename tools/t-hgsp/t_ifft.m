function [flatten_D,D] = t_ifft(flatten_A_fft,shape_ten)
    M=length(shape_ten);
    num_slices = prod(shape_ten(3:end)); % n3 x n4 x ... x np
    flatten_shape_D = [shape_ten(1:2), num_slices];
    
    % Computes the ifft on diffeent directions
    D = reshape(flatten_A_fft,shape_ten);
    for i=M:-1:3 %in reversed(range(2, p)):
        D = ifft(D,[], i);
    end
    flatten_D = reshape(D,flatten_shape_D);

    
    
    
    
    
  

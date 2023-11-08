function [SpatialOut,FourierOut] = t_eig(A,ordereig)
    % Take a tensors of order p, output their t eigendecomposition. 
    % Suppose A is in dimension n1 x n2 x ... x np.
    % ordereig = Largest Magnitude 'LM' or Smallest Magnitude 'SM' 
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
    end
    
    flatten_V_fft = zeros(flatten_shape_D);
    flatten_S_fft = zeros(flatten_shape_D);
    flatten_Yt = zeros(flatten_shape_D);
    flatten_Yt2 = zeros(flatten_shape_D);

    %%% Regular eigen value decomposition for each slices
    for i=1:num_slices
        if ~issymmetric(flatten_D(:,:,i))
            pause
        end
        rng(0)
        [flatten_V_fft(:,:,i),flatten_S_fft(:,:,i)]= eigs(flatten_D(:,:,i),flatten_shape_D(1),ordereig);
        %[flatten_V_fft(:,:,i),flatten_S_fft(:,:,i)]= eig(flatten_D(:,:,i));
        maxEig=max(max(abs(flatten_S_fft(:,:,i))));
        flatten_Yt(:,:,i) = flatten_V_fft(:,:,i)* abs(flatten_S_fft(:,:,i))/maxEig *flatten_V_fft(:,:,i)';
        flatten_Yt2(:,:,i) = flatten_V_fft(:,:,i)* (flatten_S_fft(:,:,i)/maxEig) * flatten_V_fft(:,:,i)';
    end
    %%% Unfold all dimensions after 2
    S = reshape(flatten_S_fft,shape_ten);
    V = reshape(flatten_V_fft,shape_ten);
    
    
    %     powerIm=np.sum(S.imag**2)
    %     if powerIm<1e-10:
    %         S=S.real
    %     powerIm=np.sum(V.imag**2)
    %     if powerIm<1e-10:
    %         V=V.real
    %     

    Shift_Operator_Abs_Norm = reshape(flatten_Yt,shape_ten);
    Shift_Operator_Norm = reshape(flatten_Yt2,shape_ten);
    
    
    for i=p:-1:3 %in reversed(range(2, p)):
        S = ifft(S,[], i,'symmetric');
        V = ifft(V,[], i,'symmetric');
        Shift_Operator_Abs_Norm = ifft(Shift_Operator_Abs_Norm,[], i);
        Shift_Operator_Norm = ifft(Shift_Operator_Norm,[], i);
    end
    powerIm=sum(imag(S(:)).^2);
    if powerIm<tol
        S=real(S);
    end
    powerIm=sum(imag(V(:)).^2);
    if powerIm<tol
        V=real(V);
    end
    
    SpatialOut.V=V;
    SpatialOut.S=S;
    SpatialOut.Shift_Operator_Abs_Norm=Shift_Operator_Abs_Norm;
    SpatialOut.Shift_Operator_Norm=Shift_Operator_Norm;
    
    FourierOut.V=flatten_V_fft;
    FourierOut.S=flatten_S_fft;
    
    
    
    
    
  

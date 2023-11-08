function H = HyperRepresentation(B)
    % B is the incidence matrix
    N=size(B,1);
    M=max(sum(B));
    A = adjacency_tensor(B,N,M); 
    
    % Degree tensor
    d=zeros(N,1);
    for i=1:N
        otherdims = repmat({':'},1,ndims(A)-1);
        d(i) = sum(double(A(i,otherdims{:})),'all');
    end
    X.size=size(A);
    X.subs=repmat(1:N,M,1)';
    X.vals=d;
    D = sptensor(X.subs,X.vals,X.size) ;
    L=D-A;
    H.A=A;
    H.B=B;
    H.D=D;
    H.L=L;
    
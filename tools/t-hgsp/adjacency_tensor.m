function A = adjacency_tensor(H,N,M)
    sizeTen = repmat(N,1,M);
    A=sptensor([],[],sizeTen);
    %generate the Mth order N dimension adjacency tensor from the incidence matrix H.
    for i=1:size(H,2)
        location = find(H(:,i) ~= 0);
        all_perms = generate_perms(location,M);  % all possible permutations of indices

        num_of_all_perms = size(all_perms,1); 
        c = length(location);
        A(all_perms) = c / num_of_all_perms;
    end
end
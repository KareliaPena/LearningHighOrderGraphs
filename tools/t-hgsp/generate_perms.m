function final_lst = generate_perms(hyperedge,M)
    %Given the list of nodes contrained in a hyperedge, and M, generate all
    %possible permutations of the indices.
    if length(hyperedge) == M 
        
        final_lst=perms(hyperedge);
        
    elseif length(hyperedge) == 2 && M==3
        complementary_lst_part1 = hyperedge; % the M-c indices
        complementary_lst=[repmat(hyperedge',size(complementary_lst_part1,1),1),complementary_lst_part1];
        final_lst=[];
        for j=1:size(complementary_lst,1)
            C = perms(complementary_lst(j,:));
            final_lst=[final_lst;C];
        end
        final_lst = unique(final_lst , 'rows');
        
    else
        %%Only consider the case that c < M.
        M_minus_c = M - length(hyperedge);
        complementary_lst_part1 = combinations_with_replacement(hyperedge, M_minus_c); % the M-c indices
        complementary_lst=[repmat(hyperedge',size(complementary_lst_part1,1),1),complementary_lst_part1];
        final_lst=[];
        for j=1:size(complementary_lst,1)
            C = perms(complementary_lst(j,:));
            final_lst=[final_lst;C];
        end
        final_lst = unique(final_lst , 'rows');
    end
end

function C = combinations_with_replacement(lst, M_minus_c)
    [args{1:M_minus_c}] = deal(lst);
    n = length(args);

    [F{1:n}] = ndgrid(args{:});

    for i=n:-1:1
        G(:,i) = F{i}(:);
    end

    C = unique(G , 'rows');
end
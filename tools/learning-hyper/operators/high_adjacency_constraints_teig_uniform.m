function [mat_obj_ifft_tube_dir,R,P,K,norm_R,normPtP,num_of_all_perms,listUnique,NDE]=high_adjacency_constraints_teig_uniform(N,M)
%% matrix for objective (vech -> vec) 
% takes nonsymmtric elements of the cube organized in frontal slices 
% and returns a complete cube of (N+1) frontal slices

listUnique = [];
v=1:N;
NDE=0;
for i=M
    NDE=NDE+nchoosek(N,i);
    C = num2cell(nchoosek(v,i),2);
    listUnique=[listUnique;C];
end

disp('Step 1/7')


%% Operator that accounts for symmetric: from vector to whole tensor
tic
dims = repmat(N,1,M);
cell_math=cell(NDE,1);
parfor ii=1:NDE
    elm=listUnique{ii};
    all_perms = generate_perms(elm',M);  % all possible permutations of indices
    num_of_all_perms(ii) = size(all_perms,1);
    all_perms_cell=cell(M,1);
    for j=1:M
        all_perms_cell{j}=all_perms(:,j);
    end
    ind = sub2ind(dims,all_perms_cell{:});
    cell_math{ii,1}=[ind,repmat(ii,length(ind),1)];
end
max_alpha=max(num_of_all_perms);
elem_sparse=cell2mat(cell_math);
mat_obj=sparse(elem_sparse(:,1),elem_sparse(:,2),ones(size(elem_sparse,1),1),N^M,NDE);
clear elem_sparse cell_math
disp('Step 2/11')
time=toc
%% Degree constraint

temp=ones(1,N^(M-1));
ACell = repmat({temp}, 1, N);
mat_obj_degree = sparse(blkdiag(ACell{:}));

clear ACell
disp('Step 3/11')

%% This math object accounts for the symmetry in the third dimenssion 

mat_obj2_cell=cell(1,1);
shape_As=[N,N,repmat(2*N+1,1,M-2)];
for m=3:M
    frontal=sparse(prod(shape_As(1:m-1)),prod(shape_As(1:m-1))*N);
    top=speye(prod(shape_As(1:m-1))*N);
    ACell_1 = repmat({sparse(fliplr(speye(prod(shape_As(1:m-1)))))}, 1, N);
    bottom = fliplr(matlab.internal.math.blkdiag(ACell_1{:}));
    mat_obj2_p1=1/2*sparse([frontal;top;bottom]);
    if m<M
        ACell_2 = repmat({mat_obj2_p1}, 1, N^(M-m));
        mat_obj2_p1 = matlab.internal.math.blkdiag(ACell_2{:});
    end
    mat_obj2_cell{m-2} = mat_obj2_p1;
end

mat_obj2=mat_obj2_cell{1};
for j=2:length(mat_obj2_cell)
    mat_obj2=mat_obj2_cell{j}*mat_obj2;
end

 
clear frontal top bottom mat_obj2_cell ACell_2 mat_obj2_p1 ACell_1

disp('Step 4/11')


%% Everything together to get the fft on the different directions
mat_obj_fft = real(dftmtx(2*N+1));
ACell = repmat({mat_obj_fft}, 1, prod(shape_As(1:M-1)));
mat_obj_fft_block =matlab.internal.math.blkdiag(ACell{:});


mat_obj_fft_all_dir=1;
for ii=3:M
    i=1:prod(shape_As(ii:M));
    j=[];
    for p=1:prod(shape_As(1:ii-1))
        j=[j,(i-1)*prod(shape_As(1:ii-1))+p]; % select here the step 
    end
    mat_obj3=sparse(1:prod(shape_As),j,ones(length(j),1),prod(shape_As),prod(shape_As));
    %tSymTensor_tubal{ii-2}=reshape( mat_obj3{ii-2}*As_full(:),2*N+1,[]);
    mat_obj_fft_all_dir=mat_obj3'*mat_obj_fft_block*mat_obj3*mat_obj_fft_all_dir;
end


disp('Step 5/11')


%% Everything together to get the ifft on the different directions
shape_As_tube=[1,1,shape_As(3:end)];
mat_obj_ifft = conj(mat_obj_fft)/(2*N+1);
ACell = repmat({mat_obj_ifft}, 1, prod(shape_As_tube(1:M-1)));
mat_obj_ifft_block =matlab.internal.math.blkdiag(ACell{:});

mat_obj_ifft_tube_dir=1;
for ii=3:M
    i=1:prod(shape_As_tube(ii:M));
    j=[];
    for p=1:prod(shape_As_tube(1:ii-1))
        j=[j,(i-1)*prod(shape_As_tube(1:ii-1))+p]; % select here the step 
    end
    mat_obj4=sparse(1:prod(shape_As_tube),j,ones(length(j),1),prod(shape_As_tube),prod(shape_As_tube));

    %tSymTensor_tubal{ii-2}=reshape( mat_obj4{ii-2}*As_tube(:),2*N+1,[]);
    mat_obj_ifft_tube_dir=mat_obj4'*mat_obj_ifft_block*mat_obj4*mat_obj_ifft_tube_dir;
end

disp('Step 6/11')

%% Only the operators that I need
tic
R = mat_obj_degree*mat_obj;
P = mat_obj2*mat_obj; %-> M_dup3s. Mdup
mat_obj_ifft = conj(mat_obj_fft)/(2*N+1);
K =mat_obj_fft_all_dir; %-> M_t2f. T. M_f2t
toc
disp('Step 7/11')
%% Compute the norm of R
%TODO: estimate this norm in terms of M and N
tic
norm_R = normR(R);
normPtP=max_alpha/2^(M-2);
toc
disp('Step 8/11')

%%
folder='G:\Shared drives\Learning_Hypergraphs\Operators';
save([folder '\A_uniform_M=' num2str(M) 'N=' num2str(N)],'normPtP','norm_R','mat_obj_ifft_tube_dir','R','P','mat_obj_ifft','K','num_of_all_perms','NDE','listUnique','-v7.3');
end


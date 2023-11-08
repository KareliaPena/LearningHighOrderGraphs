function laplacian_constraints_teig_spatial_uniform(N,M)
%% matrix for objective (vech -> vec) 
% takes nonsymmtric elements of the cube organized in frontal slices 
% and returns a complete cube of (N+1) frontal slices

disp('Computing operators')

listUnique = [];
v=1:N;
NDE=0;
for i=[1 M]
    NDE=NDE+nchoosek(N,i);
    C = num2cell(nchoosek(v,i),2);
    listUnique=[listUnique;C];
end

equal_list=1:N

% NDE=nchoosek(N,M)+nchoosek(N,2)+N;
% list = cell(NDE,1);
% listUnique = cell(NDE,1);
% kk=0;
% equal_list=[];
% for i=1:N
%     for j=1:N
%         for k=1:N
%             if i==j && j==k
%                 kk=1+kk;
%                 list{kk}=[i,j,k];
%                 C = unique(list{kk});
%                 listUnique{kk}=C;
%                 equal_list=[equal_list,kk];
%             else
%                 if i<=j && j<k
%                     kk=1+kk;
%                     list{kk}=[i,j,k];
%                     C = unique(list{kk});
%                     listUnique{kk}=C;
%                 end    
%             end
%             
%         end
%     end
% end
disp('Step 1/11')


%% Operator that accounts for symmetric: from vector to whole tensor
%MATH_Cell=cell(NDE,1);
mat_obj=sparse(N^M,NDE);

for ii=1:NDE
    elm=listUnique{ii};
    if length(elm)==1
        C=[elm,elm,elm];
    elseif length(elm)==2
        Element1=perms([elm(1),elm(1),elm(2)]);
        Element2=perms([elm(1),elm(2),elm(2)]);
        Element=[Element1;Element2];
        [C,ia,ic] = unique(Element,'rows');
    elseif length(elm)==3
        Element=perms(elm);
        [C,ia,ic] = unique(Element,'rows');
    end
    ind = sub2ind([N,N,N],C(:,1),C(:,2),C(:,3));
%    MATH_Cell{ii,1}=sparse(1,ind+N^2,1,1,N^M+N^2) ;
    mat_obj(ind,ii)=1;
end


%mat_obj = sparse(cell2mat(MATH_Cell)');


disp('Step 2/11')


%% This math object accounts for the symmetry in the third dimenssion 
frontal=sparse(N^2,N^3);
top=speye(N^3);
bottom=sparse(N^3,N^3);
for i=1:N
    bottom(1+(i-1)*N^2:i*N^2,1+(N-i)*N^2:(N-i+1)*N^2)=speye(N^2);
end
mat_obj2=1/2*sparse([frontal;top;bottom]);

disp('Step 3/11')

% SymmetricTensor=mat_obj2*CompleteTensorVec;

% CompleteSymTensor=tensor(reshape(SymmetricTensor,[N,N,2*N+1]));

%% This math object maps the data cube organized from frontal slices to scalar tubes
i=1:2*N+1;
j=(i-1)*N^2+1;
aux=sparse(i,j,ones(length(i),1),2*N+1,N^2*(2*N+1));

mat_obj3=aux;
for i=1:N^2-1
    mat_obj3=[mat_obj3;circshift(aux,i,2)];
end

disp('Step 4/11')

% TubesTensor=mat_obj3*SymmetricTensor;

%% mat_obj_fft maps the data cube organized in scalar tubes from the spatial domain to the fft domain

aux1 = real(dftmtx(2*N+1));
mat_obj_fft=aux1;
% Elem= N^2;
% ACell = repmat({mat_obj_fft}, 1, Elem);
% mat_obj_fft = sparse(blkdiag(ACell{:}));

% ACell = repmat({mat_obj_fft}, 1, N);
% mat_obj_fft = sparse(blkdiag(ACell{:}));

% Tensorffttubes=mat_obj_fft*TubesTensor;
disp('Step 5/11')


%% This math object maps the data cube organized from scalar tubes to frontal slices
i=1:N^2;
j=(i-1)*(2*N+1)+1;

aux=sparse(i,j,ones(length(i),1),N^2,N^2*(2*N+1));

mat_obj4=aux;
for i=1:(2*N)
    mat_obj4=[mat_obj4;circshift(aux,i,2)];
end
mat_obj4=sparse(mat_obj4);


disp('Step 6/11')

% Tensorfftfrontal=mat_obj4*Tensorffttubes;





%% Math constraint that makes non super diagonal elements negative
mat_const_ndiag=speye(NDE,NDE);
for ii=1:N
    mat_const_ndiag(equal_list(ii),equal_list(ii))=0;
end

% NonDiagElemtns=mat_const_diag*NNE_test;

%% Math constraint that makes the super diagonal equal to 1
mat_const_diag=sparse(N,NDE);
for ii=1:N
    mat_const_diag(ii,equal_list(ii))=1;
end

% DiagElemtns=mat_const_diag*NNE_test;



%% Math constraint that makes sum of the frontal slices equal to zero
%% THIS MATH CONSTRAIN MIGHT NOT BE NECESSARY
temp=ones(1,N^2);
ACell = repmat({temp}, 1, N);
mat_const_sum_front = sparse(blkdiag(ACell{:}));
mat_const_sum_front_F=mat_const_sum_front*mat_obj;
% SumFront=mat_const_sum_front*CompleteTensorVec;



%% create constraint matrices
% equality constraint A2*vech(L)==b1
A1 = [mat_const_sum_front_F;mat_const_diag];
b1 = [sparse(N,1);ones(N,1)];


% inequality constraint A1*vech(L)<=b2
A2 = mat_const_ndiag;
b2 = sparse(NDE,1);

disp('Step 11/11')
folder='G:\Shared drives\Learning_Hypergraphs\Operators\';
save([folder '\L_uniform_N=' num2str(N)],'equal_list','A1','b1','A2','b2','mat_obj','mat_obj2','mat_obj3','mat_obj4','mat_obj_fft','NDE','listUnique','-v7');
end
function [vec_As, vech_As, listUnique,stat] = high_uniform_pdl_hgsp(Zs_fft, alpha, beta ,params)
%high_uniform_pdl_hgsp Learn hypergraph from pairwise distances
%
%   Usage:  [vec_As, stat] = high_uniform_pdl_hgsp(Zs_fft, alpha, beta)
%           [vec_As, stat] = high_uniform_pdl_hgsp(Zs_fft, alpha, beta ,params)
%
%   Inputs:
%         Zs_fft    : Pair-wise distance tensor
%         alpha     : Log prior constant  (bigger alpha -> bigger weights in v)
%         beta      : Controls the sparsity
%         params    : Optional parameters
%
%   Outputs:
%         vec_As     : Vector containing all the elements of the adjacency tensor 
%         stat       : Optional output statistics (adds small overhead)
%
%   'vec_As = high_uniform_pdl_hgsp(Z, a, params)' computes the vector form 
%   of the adjacency tensor from the squared distances between pairs in Z, using the 
%   smoothness assumption that X^TLX is small, where X is the data (columns) 
%   that change smoothly from node to node in the hypergraph.
%
%
%   The minimization problem solved is 
%
%      argmin_vech(As) aggregate(combine(As * Zs)) - alpha * sum(log( R * vech_As))) + beta * ||As||_F^2
%
%   Where:
%      aggregate(combine(As * Zs)) and ||As||_F^2 can can be expressed 
%      as a function of vech_As (see references).
%
%   The algorithm used is forward-backward-forward (FBF) based primal dual
%   optimization (see references).
%
%%   Example:
% 
%   Size = 3; M = 3;
%   I = [1,0,1; 0,1,0; 1,0,1];
%   [coordsX,coordsY] = meshgrid(1:size(I,2),1:size(I,1));
%   coords = [coordsX(:),coordsY(:)];
%   Gray = I(:);
%   N = length(Gray);
% 
%   #Build hypergraph signal
%   [X_tensor] = tensorSignal(Gray,M);
%   Xs = double(symmetrize_tensor(X_tensor));
%   Xs_fft = fft(Xs,[],3);
%   powerIm = sum(sum(sum(imag(Xs_fft).^2)));
%   if powerIm < eps
%     Xs_fft = real(Xs_fft);
%   end
% 
%   #Build distance hypergraph signal Compute distances for each frontal slice independently
%   Zs_fft=zeros(N,N,2*N+1);
%   for k=1:size(Xs_fft,3)
%       Zs_fft(:,:,k) = gsp_distanz(Xs_fft(:,:,k)').^2;
%   end
%   params.verbosity = 5;
%
%   #Hypergraph Learning
%   alpha = 1; beta = 1;
%   [vec_As, stat] = high_uniform_pdl_hgsp(Zs_fft, alpha, beta,params);
%       
%
%%   Additional parameters
%   ---------------------
%
%    verbosity       : Above 1 adds a small overhead. Default: 1
%    maxit           : Maximum number of iterations. Default: 1000
%    tol             : Tolerance for stopping criterion. Default: 1e-5
%    step_size       : Step size from the interval (0,1). Default: 0.5
%    max_v           : Maximum weight allowed for each edge (or inf)
%
%   The stopping criterion is whether both relative primal and dual
%   distance between two iterations are below a given tolerance. 
%   
%   To set the step size use the following rule of thumb: Set it so that
%   relative change of primal and dual converge with similar rates (use
%   verbosity > 1).
%
%
%   See also: gsp_learn_graph_log_degrees gsp_distanz gsp_update_weights
% 
%   References:
%     V. Kalofolias. How to learn a graph from smooth signals. Technical
%     report, AISTATS 2016: proceedings at Journal of Machine Learning
%     Research (JMLR)., 2016.
%     
%     N. Komodakis and J.-C. Pesquet. Playing with duality: An overview of
%     recent primal? dual approaches for solving large-scale optimization
%     problems. Signal Processing Magazine, IEEE, 32(6):31--54, 2015.
%     
%     V. Kalofolias and N. Perraudin. Large Scale Graph Learning from Smooth
%     Signals. arXiv preprint arXiv:1710.05654, 2017.
%
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
%
%     K. Pena-pena, L. Taipe, F. Wang: Learning Hypergraphs from Data 
%     Using t-Products. 14(8), 1–11, (2021).
%
%
% Author: Karelia Pena and Lucas Taipe 
% Testing: *****************
% Date: December 21


if nargin < 3
    params = struct;
end

%% Default parameters
if not(isfield(params, 'verbosity')),   params.verbosity = 1;   end
if not(isfield(params, 'maxit')),       params.maxit = 1000;      end
if not(isfield(params, 'tol')),         params.tol = 1e-5;      end
if not(isfield(params, 'step_size')),   params.step_size = .5;      end     % from (0, 1)
if not(isfield(params, 'max_v')),       params.max_v = inf;         end
      

%% Fix parameter size and initialize
M = length(size(Zs_fft)); % hypergraph order
N = size(Zs_fft,1); % Number of nodes

%% Operators:
% K : Operator to pass a form vector to fourier space
% P : Operator that builds a tensor from the unique values
% R : Operator to obtain the unique elements of the degree tensor from the adjacency tensor, the operation is reversible.
% mat_obj_ifft_tube_dir : Discrete Fourier transform matrix is equivalent to \Gamma


if not(isfield(params, 'operators'))
    try
        load(['A_uniform_M=' num2str(M) 'N=' num2str(N)])
        %load(['A_M=' num2str(M) 'N=' num2str(N)])
    catch
        [mat_obj_ifft_tube_dir, R, P, K, norm_R, normPtP, num_of_all_perms, listUnique, NDE] = high_adjacency_constraints_teig_uniform(N,M);
    end
    shape_Zs = size(Zs_fft);
    num_slices = prod(shape_Zs(3:end));
    Zs_fft_flattened = reshape(Zs_fft,[N,N,num_slices]);
    disp('Computing CKP')
    CKP = zeros(num_slices,NDE);
    KP = K*P;
    for i=1:num_slices
        vecZsk = Zs_fft_flattened(:,:,i);
        CKP(i,:) = vecZsk(:)'*KP(N^2*(i-1)+1:N^2*i,:);
    end
else
    R = params.operators.R;
    CKP = params.operators.CKP; % the operator CKP contains of the operator J
    P = params.operators.P;
    norm_R = params.operators.norm_R;
    mat_obj_ifft_tube_dir = params.operators.mat_obj_ifft_tube_dir;
    listUnique = params.operators.listUnique;
end

% NDE number of unique elements
NDE = length(listUnique);

%Vector form que contine los calores unicos del tensor A_s
vech_As = zeros(NDE,1);

%% Learn the hypergraph
% min_vech_As    f(vech_As)          +       g(R_op(vech_As))      +   h(vech_As
% min_vech_As    sum(mat_obj_ifft_tube_dir*CKP*vech_As)  - alpha * sum(log(R_op(vech_As))) +beta*(P*vech_As)'*(P*vech_As)

%To obtain FBF(for our model), we need the following:
f.eval = @(vech_As) sum(mat_obj_ifft_tube_dir*CKP*vech_As);
f.prox = @(vech_As, c) min(params.max_v, max(0, vech_As - c*(ones(1,size(mat_obj_ifft_tube_dir,1))*mat_obj_ifft_tube_dir*CKP)'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_prox_log.verbose = params.verbosity - 3;
g.eval = @(z) -alpha * sum(log(z));
g_star_prox = @(z, c) z - c*alpha * prox_sum_log(z/(c*alpha), 1/(c*alpha), param_prox_log);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h.eval = @(vech_As) beta*(P*vech_As)'*(P*vech_As);
h.grad = @(vech_As) 2*beta*(P'*P)*vech_As;
% Lipschitz constant of h.grad:
h.beta = 2 * beta * normPtP;
% R: edges -> nodes
R_op = @(v) R*v;
% R': nodes -> edges
Rt_op = @(z) R'*z;


%% Custom FBF based primal dual (see [1] = [Komodakis, Pesquet])
% parameters mu, epsilon for convergence (see [1])
mu = h.beta + norm_R;     
epsilon = lin_map(0.0, [0, 1/(1+mu)], [0,1]);   % in (0, 1/(1+mu) )
    
% INITIALIZATION
% primal variable ALREADY INITIALIZED
% v : primal variable
% d : dual variable

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aprroach in hypergraph
v = vech_As;

d_n = R_op(v);
if nargout > 1 || params.verbosity > 1
    stat.f_eval = nan(params.maxit, 1);
    stat.g_eval = nan(params.maxit, 1);
    stat.h_eval = nan(params.maxit, 1);
    stat.fgh_eval = nan(params.maxit, 1);
    stat.pos_violation = nan(params.maxit, 1);
end
if params.verbosity > 1
    fprintf('Relative change of primal, dual variables, and objective fun\n');
end
tic
gn = lin_map(params.step_size, [epsilon, (1-epsilon)/mu], [0,1]); % in [epsilon, (1-epsilon)/mu]

for i = 1:params.maxit
    
    Y_n = v - gn * (h.grad(v) + 2*Rt_op(d_n));
    y_n = d_n + gn * (R_op(v));
    P_n = f.prox(Y_n, gn);
    p_n = g_star_prox(y_n, gn);
    Q_n = P_n - gn * (h.grad(P_n) + 2*Rt_op(p_n));
    q_n = p_n + gn * (R_op(P_n));
    
    if nargout > 1 || params.verbosity > 2
        stat.f_eval(i) = f.eval(v);
        stat.g_eval(i) = g.eval(R_op(v));
        stat.h_eval(i) = h.eval(v);
        stat.fgh_eval(i) = stat.f_eval(i) + stat.g_eval(i) + stat.h_eval(i);
        stat.pos_violation(i) = -sum(min(0,v));
        stat.degrees{i} = R_op(v);
    end
    
    rel_norm_primal = norm(- Y_n + Q_n, 'fro')/norm(v, 'fro');
    rel_norm_dual = norm(- y_n + q_n)/norm(d_n);
    
    if params.verbosity > 2
        fprintf('iter %4d: %6.4e   %6.4e   %6.4e    %6.4e   %6.3e \n', i, rel_norm_primal, rel_norm_dual, stat.fgh_eval(i), full(sum(v)*2), stat.h_eval(i));
    elseif params.verbosity > 1
        fprintf('iter %4d: %6.4e   %6.4e\n', i, rel_norm_primal, rel_norm_dual);
    end
    v = v - Y_n + Q_n;
    d_n = d_n - y_n + q_n;
    
    if rel_norm_primal < params.tol && rel_norm_dual < params.tol
        break
    end
end
stat.time = toc;
if params.verbosity > 0
    fprintf('# iters: %4d. Rel primal: %6.4e Rel dual: %6.4e   %6.3e\n', i, rel_norm_primal, rel_norm_dual, f.eval(v) + g.eval(R_op(v)) + h.eval(v));
    fprintf('Time needed is %f seconds\n', stat.time);
end

if params.verbosity > 3
    figure; plot(real([stat.f_eval, stat.g_eval, stat.h_eval])); hold all; plot(real(stat.fgh_eval), '.'); legend('f', 'g', 'h', 'f+g+h');
    figure; plot(stat.pos_violation); title('sum of negative (invalid) values per iteration')
    figure; semilogy(max(0,-diff(real(stat.fgh_eval'))),'b.-'); hold on; semilogy(max(0,diff(real(stat.fgh_eval'))),'ro-'); title('|f(i)-f(i-1)|'); legend('going down','going up');
end
% We obtain the adjacency tensor in its vector form
vech_As = v;
vec_As = P*v;
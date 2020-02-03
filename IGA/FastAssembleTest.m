clear all
clc
tic
load circplate_coarse.mat
Model.DegreeElevate(5,1);
Model.DegreeElevate(5,2);
Model.DegreeElevate(6,3);
Model.KnotRefine(1/3:1/3:1-1/3,1);
Model.KnotRefine(1/3:1-1/3:1-1/3,2);
Model.KnotRefine(1/3:1-1/3:1-1/3,3);
p = [Model.pu; Model.pv; Model.pw];
knots = {Model.U; Model.V; Model.W};
Rules = getWQ(p,knots);
% Jacobian of quadrature points. turning warnings off
% warning('off','all');
 [Basis, DerBasis, qpoints, J, ControlPoints] = coordinates_and_jacobian(Model,Rules);
% warning('on','all');
% Evaluate c(x) on quadrature points
% For mass matrix:
nqu = length(Rules{1}.Points);
nqv = length(Rules{2}.Points);
nqw = length(Rules{3}.Points);

%% Stiffness:
YOUNG = 30*10^9;
POISSON = 0.2;
D = get_matprop_matrix(1,YOUNG, POISSON);
% C = cell(Rules{1}.ndof,Rules{2}.ndof,Rules{3}.ndof,3);
c = reshape(c,[size(qpoints),1,1]);
[q1, q2, q3] = size(qpoints);
A1 = zeros(Rules{3}.ndof,1,q1,q2);
A2 = zeros(Rules{3}.ndof*Rules{2}.ndof,1,q1);
A3 = zeros(Rules{3}.ndof*Rules{2}.ndof*Rules{1}.ndof,1);
%% Sum-Factorization
% To correctly assign each row to its degree of freedom, we must build an
% index matrix, which is exactly the INN matrix from the book.
% You can run the geopdes' code and find that the "indices" matrix is the
% same thing as the INN index.
% Ctrl+C + Ctrl+V the code for reference
% d = ndims(fun_val);
% N_dof = space.ndof;
% 
% nonzeros = prod(arrayfun(@(i)sum(Connectivity(i).num_neigh),1:d));
% cols = zeros(1,nonzeros); rows = cols; values = cols;
% ncounter = 0;
% 
% n_size = space.ndof_dir;
% indices = cell(1,d);
% [indices{:}] = ind2sub(n_size, 1:N_dof);
% indices = cell2mat(indices);  indices = reshape(indices,[N_dof d]);
% points = cell(1,d); j_act = cell(1,d); len_j_act = zeros(1,d); n_index = zeros(1,d);
% for ll = 1:d
%     n_index(ll) = prod(n_size(1:ll-1));
% end
[INN, IEN, nel, nnz] = Model.get_connectivity; %nen = nnz
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN)));
cols = zeros(1,nnz); rows = cols; values = cols;
ncounter = 0;
[~, d] = size(INN);
n_size = [Rules{1}.ndof Rules{2}.ndof Rules{3}.ndof];
NDOF = prod(n_size);
n_idx = zeros(1,d);
for ll=1:d
    n_idx(ll) = prod(n_size(1:ll-1));
end

for ii=1:NDOF % Row-loop
    idx = INN(ii,:);
    for ll=1:d
        j_act{ll} = Rules{ll}.I{ind(ll)}; % Active trial functions
        len_j_act(ll) = length(j_act{ll});
    end
    i_nnz = prod(len_j_act);
    % Stage 1
    for k=1:q3
        B = Rules{3}.Basis(:,k)*Rules{3}.Basis(1,k)';
        for j=1:q2
            for i=1:q1
                A1(:,:,i,j) = A1(:,:,i,j) + c(i,j,k,1,1)*B;
            end
        end
    end
    % Stage 2
    for j=1:q2
        B = Rules{2}.Basis(:,j)*Rules{2}.Basis(1,j)';
        for i=1:q1
            A2(:,:,i) = A2(:,:,i) +kron(A1(:,:,i,j),B);
        end
    end
    % Stage 3
    for i=1:q1
        B = Rules{1}.Basis(:,i)*Rules{1}.Basis(1,i)';
        A3(:,:) = A3(:,:) +kron(A2(:,:,i),B);
    end
    values(ncounter+1:ncounter+i_nnz) = A3(:)';
    
    rows(ncounter+1:ncounter+i_nnz) = ii;
    i_col = zeros(d,i_nnz);
    for ll=1:d
        rep = len_j_act; rep(ll)=1;
        perm = ones(1,d); perm(ll) = len_j_act(ll);
        ap = repmat(reshape(j_act{ll}',perm),rep);
        i_col(ll,:) = ap(:)';
    end
    cols(ncounter+1:ncounter+i_nnz) = 1+n_idx*(i_col-1);
    ncounter = ncounter+i_nnz;
                
end
M = sparse(rows,cols,values, 3*NDOF, 3*NDOF);
toc
                

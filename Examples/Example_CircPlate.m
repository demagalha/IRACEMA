clear all
close all
clc
%% MIT License
% 
% Copyright (c) 2019 Guilherme Henrique da Silva and André Demetrio de Magalhães
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%% Circular Plate Modal Analysis
% This example is the modal analysis of a circular plate p-refined
% presented by Cottrell, Austin, Reali, Bazilev and Hughes in the
% 2006 "Isogeometric Analysis of Structural Vibrations" paper

%% Pre-processing
load circplate_coarse.mat
% Elevate pu to 4, pv to 5 and pw to 2
% u = 1st, v = 2nd, w = 3rd
Model.DegreeElevate(2,1); % augment 2 degrees in 1st parametric direction 
Model.DegreeElevate(3,2); % augment 3 degrees in 2nd parametric direction
Model.DegreeElevate(1,3); % augment 1 degree in 3rd parametric direction

% Note that our disc is parametrized differently than the article's. We use
% only one element to describe the whole geometry, while the article uses
% an eight element mesh. So we need to refine the mesh a little bit. We
% will add 2 knots in u and v directions
Model.KnotRefine([0.33 0.66],1);
Model.KnotRefine([0.33 0.66],2);
% This way, there is 9 elements in total, not 8, since our parametrization
% is slightly different. As you will notice, the autovalues will also be
% different.
% Material Properties
YOUNG_MODULUS = 30*10^9;
RHO = 2.32*10^3; 
POISSON = 0.2;
D = get_matprop_matrix(1,YOUNG_MODULUS,POISSON); % Isotropic Deformation

% Generate Connectivity arrays
[INN, IEN, nel, nen] = Model.get_connectivity;
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN)));
LM = zeros(3*nen,nel);
for i = 1 : nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),3*nen,1);
end
% Get Model parameters
pu = Model.pu;
pv = Model.pv;
pw = Model.pw;
U = Model.U;
V = Model.V;
W = Model.W;
P = Model.get_point_cell;

% Get Gauss-Legendre Quadrature Points
[u, wu] = getGP(pu);
[v, wv] = getGP(pv);
[w, ww] = getGP(pw);
N_QUAD_U = length(u);
N_QUAD_V = length(v);
N_QUAD_W = length(w);

% Alocate memory for Stiffness and Mass arrays
K = zeros(numel(INN),numel(INN));
M = K;
N_DOF = numel(INN);
N_ELE_DOF = nen*3;

%% Assembly
for e=1:nel % Loop Through Elements
    ni = INN(IEN(1,e),1);
    nj = INN(IEN(1,e),2);
    nk = INN(IEN(1,e),3);
    % Check if element has zero measure
    if (U(ni+1) == U(ni) || (V(nj+1) == V(nj)) || (W(nk+1) == W(nk)))
        continue
    end
    K_e = zeros(3*nen,3*nen);
    M_e = K_e;
    F_e = zeros(3*nen,1);
    
    for i=1:N_QUAD_U % Loop through U quadrature points
        for j=1:N_QUAD_V % Loop through V quadrature points
            for k=1:N_QUAD_W % Loop through W quadrature points
                [R, dR, J] = Shape(u(i),v(j),w(k),e,pu,pv,pw,P,U,V,W,INN,IEN);
                Jmod = J*wu(i)*wv(j)*ww(k);
                K_e = K_e + BuildKLocal(dR,Jmod,D);
                M_e = M_e + BuildMLocal(R,Jmod,RHO);
            end
        end
    end
    % ASSEMBLY
    idx = LM(:,e);
    for i=1:N_ELE_DOF
        ii = idx(i);
        for j=1:N_ELE_DOF
            jj = idx(j);
            K(ii,jj) = K(ii,jj)+K_e(i,j);
            M(ii,jj) = M(ii,jj)+M_e(i,j);
        end
    end
end

%% Apply boundary conditions
% Since B-Splines are always inside the convex hull of their control
% points, to apply a boundary condition where the outermost edge of a
% circular plate is clamped, we apply the radius formulae of the circle and
% grab all elements that are outside this radius.
constNod = [];
for i = 1:numel(P)
        if (((P{i}(1)-2)^2 + (P{i}(2)-2)^2)-4 >= 0)
            P{i}(1)
            P{i}(2)
            constNod=[constNod i];
    end
end
boundaries = reshape(ID(:,constNod),numel(ID(:,constNod)),1);
% To apply the clamped boundary condition in vibrations, we usually use the
% replacement method, removing those DOFs in the matrix and replacing for
% its prescribed value in the autovector answer. However, constantly
% tinkering with matrix sizes (specially sparses) in MATLAB drains quite a
% bit of computational power, so we will apply the BCs by simply setting
% their stiffness to "infinity" so they won't move on the answer.
% For an example on the replacement method, please look the 3D Beam example
% in the Examples folder.
for i=1:numel(boundaries)
    K(boundaries(i),boundaries(i)) = 1e30;
end

%% Post-Processing
% We transform K and M into sparse matrices only now, because MATLAB is not
% very good for assembly algorithms using sparse matrices. It is faster
% when using the full form. However, since we won't touch the matrix
% anymore, we can turn them into sparses again to use the eigs() function.
K = sparse(K);
M = sparse(M);
[autovector,ome] = eigs(K,M,4,'sm');
omega = sqrt(ome);
clearvars -except K M autovector freq ome Model ID
 save('circplate.mat')
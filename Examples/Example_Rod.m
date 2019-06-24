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
clear all
close all
clc

%% ISOGEOMETRIC ANALYSIS VIBRATION OF A ROD

%% Load Data / Create Model
P = {[0 0 0 1],[1 0 0 1]}; % Control Points for 1D Rod
Model = Geometry('curve',1,[0 0 1 1],P); % 1D Rod Constructor
% Refinement Phase: First Degree Elevate, then h refine
Model.DegreeElevate(1,1); % Quadratic Elements
Model.KnotRefine([0.01:0.01:0.99],1); % Add 9 Knot Spans
P = Model.get_point_cell;
U = Model.U;
% Connectivity Arrays
[INN, IEN, nel, nen] = Model.get_connectivity;
ID = reshape(1:max(max(IEN)),1,max(max(IEN)));
LM = zeros(nen,nel);
for i = 1:nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),nen,1);
end
%% Assembly
[x,wx] = getGP(Model.pu); % Quadrature Rule (Gauss Legendre)
N_QUAD_X = length(x);
pu = Model.pu;
rho = 1; % Density
E = 1; % Elastic 
K = sparse(numel(INN),numel(INN));
M = K;
F = sparse(numel(INN),1);
N_DOF = numel(INN);
N_ELE_DOF = nen;

for e=1:nel % Loop Through Elements
    ni = INN(IEN(1,e),1);
    % Check if element has zero measure
    if U(ni+1) == U(ni)
         continue
    end
    K_e = zeros(nen,nen);
    M_e = K_e;
    F_e = zeros(nen,1);
    for i=1:N_QUAD_X % Loop through quadrature points
        [R, dR, J] = Shape1D(Model,x(i),e,P,IEN,INN);
        Jmod = J*wx(i);
        for ii=1:nen
            for jj=1:nen
                K_e(ii,jj) = K_e(ii,jj) + Jmod*dR(jj)*E*dR(ii);
                M_e(ii,jj) = M_e(ii,jj) + Jmod*R(jj)*rho*R(ii);
            end
        end
    end
    idx = LM(:,e);
    for i=1:N_ELE_DOF
        ii = idx(i);
        for j=1:N_ELE_DOF
            jj = idx(j);
            K(ii,jj) = K(ii,jj) +sparse(K_e(i,j));
            M(ii,jj) = M(ii,jj) +sparse(M_e(i,j));
        end
     end
end
%% Boundary Conditions
boundary = zeros(2,1); % Vector with boundary coordinates
boundary(1) = 0;
boundary(2) = 1; 
constNod = [];
for i=1:numel(P)
    if (P{i}(1) == boundary(1) || P{i}(1) == boundary(2))
        constNod = [constNod i];
    end
end
bc = reshape(ID(:,constNod),numel(ID(:,constNod)),1);
bc = sort(bc,'descend');
for i=1:numel(bc)
    K(:,bc(i)) = [];
    K(bc(i),:) = [];
    M(:,bc(i)) = [];
    M(bc(i),:) = [];
end
%% Post-Processing
[autovector,omega] = eigs(K,M,N_DOF-numel(bc));
bc = sort(bc,'ascend');
[sz1,sz2] = size(autovector);
boundaries = zeros(1,sz2);
bc_autos = zeros(sz1+numel(bc),sz2);
for i=1:numel(bc)
    autovector = [autovector(1:bc,:); boundaries; autovector(bc+1:end,:)];
end
omega = sqrt(diag(omega));
[omega, idx] = sort(omega,'ascend');
autovector = autovector(:,idx);
omega_t = zeros(size(omega)); % Theoretical value
for i=1:length(omega)
    omega_t(i) = i*pi;
    coordinate(i) = i/length(omega); 
end
for i=1:numel(omega)
    fun_plot(i) = omega(i)/omega_t(i);
end
plot(coordinate,fun_plot)


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
%% Assemble Function
% A simple assembly algorithm that returns the Stiffness and Mass matrices
% and the Force Array. (Force array currently not implemented)
% INPUTS ------------------------------------------------------------------
% Model             - The Geometry Object to be assembled
% MatPropMatrix     - The Materials Property Matrix. This function will be
% incorporated inside this algorithm in future versions, so it works for
% more than one material and deformation field
% (isotropic/orthotropic/anisotrophic). Currently only isotropic supported.
% OUTPUTS -----------------------------------------------------------------
% K                 - Stiffness Matrix (in sparse form)
% M                 - Mass Matrix (in sparse form)
% F                 - Force Array (currently disabled)
%--------------------------------------------------------------------------
% This function only works for Volumetric Patches (u,v,w directions).
function [K, M, IEN] = Assemble(Model,MatPropMatrix,RHO)
    D = MatPropMatrix;
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
%     check1 = min(INN);
%     check2 = max(INN);    
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
        M_e = zeros(nen,3);
%         cond1 = (ni == check1(1) || nj == check1(2) || nk == check1(3));
%         cond2 = (ni == check2(1) || nj == check2(2) || nk == check2(3));
%         if cond1 || cond2
%             % If on the boundary, use Gauss-Legendre Quadrature
%             [u, wu] = getGP(pu);
%             [v, wv] = getGP(pv);
%             [w, ww] = getGP(pw);
%         else
% % %             % If on the interior, use Cauchy-Galerkin Points for Quadrature
%             [u, wu] = getCG(pu);
%             [v, wv] = getCG(pv);
%             [w, ww] = getCG(pw);
%         end

%         F_e = zeros(3*nen,1); Currently disabled.
        for i=1:N_QUAD_U % Loop through U quadrature points
            for j=1:N_QUAD_V % Loop through V quadrature points
                for k=1:N_QUAD_W % Loop through W quadrature points
                    [R, dR, J] = Shape3D(Model,u(i),v(j),w(k),e,P,IEN,INN);
                    Jmod = abs(J*wu(i)*wv(j)*ww(k));
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
    K = sparse(K);
    M = sparse(M);
end
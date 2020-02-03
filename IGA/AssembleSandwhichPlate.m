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
%% AssembleSandwhich Function
% A simple assembly algorithm that returns the Stiffness and Mass matrices
% and the Force Array of a Sandwhich Reissner-Mindlin Plate.
%(Force array currently not implemented)
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
function [K, M, IEN] = AssembleSandwhichPlate(Model,MaterialProperties,RHO)
[INN, IEN, nel, nen] = Model.get_connectivity;
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN)));
LM = zeros(3*nen, nel);
for i = 1:nel
    LM(:,i) = reshape(ID(:,IEN(:,i)),3*nen,1);
end
% Model parameters
pu = Model.pu;
pv = Model.pv;
U = Model.U;
V = Model.V;
P = Model.get_point_cell;
% Quadrature points
[u, wu] = getGP(pu);
[v, wv] = getGP(pv);
N_QUAD_U = length(u);
N_QUAD_V = length(v);

% Memory allocation for Stiffness and Mass arrays
K = zeros(numel(INN), numel(INN));
M = K;
N_DOF = numel(INN);
N_ELE_DOF = nen*SOMETHING_I_MUST_RESEARCH_AND_CHANGE; %% plate dofs

%% Assemblàge (French for Assembly)
 for e=1:nel % Loop through elements
     ni = INN(IEN(1,e),1);
     nj = INN(IEN(1,e),2);
     % Check is element's measure is zero
     if (U(ni+1) == U(ni)) || (V(nj+1) == V(nj))
         continue
     end
     K_e = zeros(SOMETHING_I_MUST_CHANGE,SOMETHING_I_MUST_CHANGE);
     M_e = K_e;
     
     for i=1:N_QUAD_U % Loop through U quad points
         for j=1:N_QUAD_V % Loop through V quad points
             [R, dR, d2R, J] = Sandwhich3D(Model,u(i),v(j),e,P,IEN,INN);
             Jmod = J*wu(i)*wv(j);
             K_e = K_e + BuildPlateStiffness(dR,d2R,Jmod,rho_vector); %%% MUST CHANGE
             M_e = M_e + BuildPlateMass(R,Jmod,rho_vector); %%% MUST CHANGE
         end
     end
     % Assembly algorithm
     idx = LM(:,e);
     for i=1:N_ELE_DOF
         ii = idx(i);
         for j=1:N_ELE_DOF
             jj = idx(j);
             K(ii,jj) = K(ii,jj) +K_e(i,j);
             M(ii,jj) = M(ii,jj) +M_e(i,j);
         end
     end
 end
 K = sparse(K);
 M = sparse(M);     

end

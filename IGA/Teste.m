clear all
close all
clc
%% Load data
tic;
f = waitbar(0,'Carregando Dados');
pause(.5)
cd ..
cd Examples
load('beam_coarse.mat');
% h refine the beam in u,v and w directions
Model.KnotRefine([[50:50:450],[550:50:950]],1);
Model.KnotRefine([[50:50:450],[550:50:950]],2);
Model.KnotRefine(500,3);
cd ..
cd IGA
% Model = viga3;
% clear teste
[INN, IEN, nel, nen] = Model.get_connectivity;
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN)));
LM = zeros(3*nen,nel); % ~= EDOF' (Similar to EDOF in CALFEM, but transposed)
for i = 1 : nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),3*nen,1);
end
% [wx, x] = GetOptimalQuadPoints(Model.U, Model.pu);
% [wy, y] = GetOptimalQuadPoints(Model.V, Model.pv);
% [wz, z] = GetOptimalQuadPoints(Model.W, Model.pw);
% wx = 2*wx;
% wy = 2*wy;
% wz = 2*wz;
% x = 2*x -1;
% y = 2*y -1;
% z = 2*z -1;
[x, wx] = getGP(Model.pu);
[y, wy] = getGP(Model.pv);
[z, wz] = getGP(Model.pw);
N_QUAD_X = length(x);
N_QUAD_Y = length(y);
N_QUAD_Z = length(z);
pu = Model.pu;
pv = Model.pv;
pw = Model.pw;
U = Model.U;
V = Model.V;
W = Model.W;
P = Model.get_point_cell;
rho = 7860; % Densidade
Y = 210*10^9; % Rigidez
vu = 0.3;
D = get_matprop_matrix(1,Y,vu); 
K = zeros(numel(INN),numel(INN));
M = K;
F = zeros(numel(INN),1);
N_DOF = numel(INN);
N_ELE_DOF = nen*3;
delete(f);
f = waitbar(1,'Carregando Dados');
pause(.5);
delete(f);
%% Loops
f = waitbar(0,'Montando Matrizes Globais');
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
    
    for i=1:N_QUAD_X % Loop through quadx
        for j=1:N_QUAD_Y % quady
            for k=1:N_QUAD_Z % quadz
                [R, dR, J] = Shape(x(i),y(j),z(k),e,pu,pv,pw,P,U,V,W,INN,IEN);
                Jmod = J*wx(i)*wy(j)*wz(k);
                K_e = K_e + BuildKLocal(dR,Jmod,D);
                M_e = M_e + BuildMLocal(R,Jmod,rho);
%                 F_e = BuildFLocal(R,dR);
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
    
    waitbar(e/nel,f,'Montando Matrizes Globais');
end
delete(f);
f = waitbar(0,'Aplicando Condições de Contorno');
loc_y_cc1 = 0;
loc_z_cc1 = 0;
loc_y_cc2 = 1;
loc_z_cc2 = 0;
constNod = [];
for i = 1:numel(P)
    cond1 = (P{i}(2) == loc_y_cc1) || (P{i}(2) == loc_y_cc2);
    cond2 = P{i}(3) == loc_z_cc2;
    if (cond1 && cond2)
        constNod=[constNod i];
    end
end
bc1 = reshape(ID(:,constNod),numel(ID(:,constNod)),1);
rolete = length(constNod)/2;
rolete = 2:3:rolete;
tmp = 1;
for i=1:numel(bc1)
%     if (tmp <= numel(rolete)) && (constNod(i) == rolete(tmp)) 
%         tmp = tmp+1;    
%         continue
%     else
      K(bc1(i),bc1(i)) = 1e30;
      waitbar(i/numel(bc1),f,'Aplicando Condições de Contorno')
%     end
end
delete(f)
K = sparse(K);
M = sparse(M);
F = sparse(F);
[autovector,ome] = eigs(K,M,20,'sm');
freq = sqrt(diag(ome))/(2*pi);
 a = toc;
%  scaling = 10;
clearvars -except K M autovector freq ome Model ID
 save('test.mat')
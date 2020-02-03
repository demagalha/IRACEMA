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
function [K, M, IEN] = FastAssembleWQ(Model)
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
    
    nu = numel(unique(U)) -1; % # of knot-spans in U
    nv = numel(unique(V)) -1; % # of knot-spans in V
    nw = numel(unique(W)) -1; % # of knot-spans in W

    P = Model.get_point_cell;

    %% Choice of quadrature points: remark 4.1

    QPu = unique(U); % Initializing Quadrature Points in u direction
    tmp = (QPu(2)-QPu(1))/pu;
    tmp2 = (QPu(end)-QPu(end-1))/pu;
    QPu2 = QPu(1)+tmp:tmp:QPu(2)-tmp; % p+1 points on first 2 knots
    QPu3 = QPu(end-1)+tmp2:tmp2:QPu(end)-tmp2; % p+1 points on last 2 knots
    QPu4 = zeros(1,numel(QPu)-3);
    for i=1:numel(QPu4)
        QPu4(i) = (QPu(i+2)-QPu(i+1))/2 +QPu(i+1); % Midpoints
    end
    QPu = [QPu, QPu2, QPu3, QPu4]'; % Concatenating knots
    QPu = sort(QPu); % Ordering in ascending order
    clearvars QPu2 QPu3 QPu4 tmp tmp2
    
    QPv = unique(V); % Initializing Quadrature Points in v direction
    tmp = (QPv(2)-QPv(1))/pv;
    tmp2 = (QPv(end)-QPv(end-1))/pv;
    QPv2 = QPv(1)+tmp:tmp:QPv(2)-tmp; % p+1 points on first 2 knots
    QPv3 = QPv(end-1)+tmp2:tmp2:QPv(end)-tmp2; % p+1 points on last 2 knots
    QPv4 = zeros(1,numel(QPv)-3);
    for i=1:numel(QPv4)
        QPv4(i) = (QPv(i+2)-QPv(i+1))/2 +QPv(i+1); % Midpoints
    end
    QPv = [QPv, QPv2, QPv3, QPv4]'; % Concatenating knots
    QPv = sort(QPv); % Ordering in ascending order
    clearvars QPv2 QPv3 QPv4 tmp tmp2
    
    QPw = unique(W); % Initializing Quadrature Points in w direction
    tmp = (QPw(2)-QPw(1))/pw;
    tmp2 = (QPw(end)-QPw(end-1))/pw;
    QPw2 = QPw(1)+tmp:tmp:QPw(2)-tmp; % p+1 points on first 2 knots
    QPw3 = QPw(end-1)+tmp2:tmp2:QPw(end)-tmp2; % p+1 points on last 2 knots
    QPw4 = zeros(1,numel(QPw)-3);
    for i=1:numel(QPw4)
        QPw4(i) = (QPw(i+2)-QPw(i+1))/2 +QPw(i+1); % Midpoints
    end
    QPw = [QPw, QPw2, QPw3, QPw4]'; % Concatenating knots
    QPw = sort(QPw); % Ordering in ascending order
    clearvars QPw2 QPw3 QPw4 tmp tmp2
    
    [QPuu, QPvv, QPww] = ndgrid(QPu, QPv, QPw); % QP grids
    QP = cell(3,1);
    QP{1} = QPuu;
    QP{2} = QPvv;
    QP{3} = QPww;
%     for i=1:6
%         plot3(QPuu(:,:,i),QPvv(:,:,i),QPww(:,:,i),'*')
%         hold on
%     end
[sz1, sz2, sz3] = size(QPu);
sz(1) = sz1;
sz(2) = sz2;
sz(3) = sz3;

    for l=1:3
        for il=1:sz(l)
            [B, dB, ~] = Shape3D(parameters);
            B0(l,il) = 
        end
        
end
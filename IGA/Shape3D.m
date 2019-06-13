function [R, dR, J] = Shape3D(Model,qu,qv,qw,element,P,IEN,INN)
%% Comments and Function I/O
% -------------------------------------------------------------------------
% Shape3D:
% Function that calculates the shape functions of a NURBS Volume
% and its derivatives.
% -------------------------------------------------------------------------
% Inputs:
% Model - An object of the Geometry class
% qu, qv, qw - Quadrature Points on parent element coordinates [-1,1]^3
% element - The element number where the calculations are being made
% P - Cell array from Model.get_point_cell
% IEN - IEN array (Cottrell, Hughes and Bazilevs' Isogeometric Analysis)
% INN - INN array (Cottrell, Hughes and Bazilevs' Isogeometric Analysis)
% IEN and INN are not necessary, since Model has get_connectivity_array
% method, but since this function is called constantly in loops, passing
% the array inside memory on the main script is more efficient.
% -------------------------------------------------------------------------
% Outputs:
% R - NURBS Basis Functions
% dR - The first derivatives of the NURBS Basis Functions
% -------------------------------------------------------------------------

%% Setup
pu = Model.pu; pv = Model.pv; pw = Model.pw;% Polynomial Orders
nen = (pu+1)*(pv+1)*(pw+1);                 % # of local basis functions
tmp = INN(IEN(1,element),:);
nu = tmp(1); nv = tmp(2); nw = tmp(3);      % NURBS coordinates / span #
clear tmp

U = Model.U; V = Model.V; W = Model.W;      % Knot Vectors

u = ((U(nu+1)-U(nu))*qu +U(nu+1)+U(nu))/2;  % Change of variables, from
v = ((V(nv+1)-V(nv))*qv +V(nv+1)+V(nv))/2;  % [-1,1] parent element space
w = ((W(nw+1)-W(nw))*qw +W(nw+1)+W(nw))/2;  % to [0,1] parameter space

R = zeros(nen,1);                           % Trivariate NURBS basis funs
dR_duvw = zeros(nen,3);                     % Trivariate NURBS derivaties
                                            % (:,1)=du; (:,2)=dv; (:,3)=dw;
Q = 0;                                      % Weight of NURBS Basis
dQ_du = 0;                                  % Weights of NURBS Derivatives
dQ_dv = 0;
dQ_dw = 0;
weight = zeros(1,nen);                      % Array of weights in column
pts = zeros(nen,3);                         % Array of pts in column order

%% B-Spline Basis and Derivatives to NURBS Basis and Derivatives
NU = DersBasisFun(nu-1,u,pu,1,U); % Basis and 1st derivative in U direction 
NV = DersBasisFun(nv-1,v,pv,1,V); % Basis and 1st derivative in V direction
NW = DersBasisFun(nw-1,w,pw,1,W); % Basis and 1st derivative in W direction
N = NU(1,:);                      % U Knot Basis
dN = NU(2,:);                     % U Knot Derivative
M = NV(1,:);                      % V Knot Basis
dM = NV(2,:);                     % V Knot Derivative
L = NW(1,:);                      % W Knot Basis
dL = NW(2,:);                     % W Knot Derivative
clear NU NV NW

%% NURBS Rational Basis Functions

% Select the points and weights that are locally supported.
% Going descend direction (nu-uu) because that's how the book does
% And going ascend direction was giving some errors.
location = 0;
for ww=0:pw
    for vv=0:pv
        for uu=0:pu
            location = location+1;
            pts(location,:) = P{nu-uu,nv-vv,nw-ww}(1:3);
            weight(location) = P{nu-uu,nv-vv,nw-ww}(4);
        end
    end
end
% Fixing order to be ascending -> flipping up down


% Calculating the sum of weights
location = 0;
for ww=0:pw
    for vv=0:pv
        for uu=0:pu
            location = location+1;
            % Total Weighting
            Q = Q + N(pu-uu+1)*M(pv-vv+1)*L(pw-ww+1)*weight(location);
            % Weight derivatives in regards to parametric coordinates
            dQ_du = dQ_du + dN(pu-uu+1)*M(pv-vv+1)*L(pw-ww+1)*weight(location);
            dQ_dv = dQ_dv + N(pu-uu+1)*dM(pv-vv+1)*L(pw-ww+1)*weight(location);
            dQ_dw = dQ_dw + N(pu-uu+1)*M(pv-vv+1)*dL(pw-ww+1)*weight(location);

        end
    end
end

% Making the rational NURBS basis
location = 0;
for ww=0:pw
    for vv=0:pv
        for uu=0:pu
            location = location+1;
            ratio = weight(location)/(Q*Q); % Rational NURBS Weight
            BASIS = N(pu-uu+1)*M(pv-vv+1)*L(pw-ww+1); % B-Spline Basis
            R(location) = BASIS*ratio*Q; % NURBS Basis
            % NURBS Basis derivatives with regards to parametric space
            dR_duvw(location,1) = dN(pu-uu+1)*M(pv-vv+1)*L(pw-ww+1) - BASIS*dQ_du;
            dR_duvw(location,2) = N(pu-uu+1)*dM(pv-vv+1)*L(pw-ww+1) - BASIS*dQ_dv;
            dR_duvw(location,3) = N(pu-uu+1)*M(pv-vv+1)*dL(pw-ww+1) - BASIS*dQ_dw;
            dR_duvw(location,:) = dR_duvw(location,:)*ratio;
        end
    end
end

%% Jacobian
 Jacobian = pts'*dR_duvw;% Jacobian Matrix. [pts']=3xnen,[dR]=nenx3 [J]=3x3
% Some mathematical definitions here:
% The definition of Jacobian is a df_i/dx_j matrix. That means each column
% Represents one physical coordinate and each row represents one
% parametric function, in the case of B-Splines. Our Jacobian vector,
% however, is flipped. Each row is representing a physical coordinate and
% each column is representing a parametric function. So, to make up for it,
% we take the transpose here.
% PS: I realize that this could've been done by the definition, as
% Jacobian = dR_duvw'*pts, but I would've missed the opportunity on making
% these notes :) Shape2D and Shape1D will come with the correct formulae.
% With some other commentaries on the determinant of rectangular matrices
% and pseudoinverses.
J_mod = det(Jacobian);
J_2 = 0.125*(U(nu+1)-U(nu))*(V(nv+1)-V(nv))*(W(nw+1)-W(nw));
J = abs(J_mod)*J_2;
dR = dR_duvw*inv(Jacobian);
end
function [R, dR, J] = Shape1D(Model,qu,element,P,IEN,INN)
%% Comments and Function I/O
% -------------------------------------------------------------------------
% Shape2D:
% Function that calculates the shape functions of a NURBS Curve
% and its derivatives.
% -------------------------------------------------------------------------
% Inputs:
% Model - An object of the Geometry class
% qu - Quadrature Points on parent element coordinates [-1,1]
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
pu = Model.pu;           % Polynomial Orders
nen = (pu+1);            % # of local basis functions
tmp = INN(IEN(1,element),:);
nu = tmp(1);       % NURBS coordinates / span #
clear tmp

U = Model.U;       % Knot Vectors

u = ((U(nu+1)-U(nu))*qu +U(nu+1)+U(nu))/2;  % Change of variables, from
                                            % [-1,1] parent element space
                                            % to [0,1] parameter space

R = zeros(nen,1);                           % Trivariate NURBS basis funs
dR_du = zeros(nen,1);                     % Trivariate NURBS derivatives
                                            % (:,1)=du; (:,2)=dv;
Q = 0;                                      % Weight of NURBS Basis
dQ_du = 0;                                  % Weights of NURBS Derivatives
weight = zeros(1,nen);                      % Array of weights in column

dx_du = 0;                                  % Derivatives of physical to
du_dx = 0;                                  % parameter coordinate/inverse
                                            % Since this is a 1D case, no
                                            % y and z coordinates required.
                                            % For a curved rod in x,y,z
                                            % coordinates, we encourage you
                                            % to modify this routine.

%% B-Spline Basis and Derivatives to NURBS Basis and Derivatives
NU = DersBasisFun(nu-1,u,pu,1,U); % Basis and 1st derivative in U direction 
N = NU(1,:);                      % U Knot Basis
dN = NU(2,:);                     % U Knot Derivative
clear NU

%% NURBS Rational Basis Functions

% Select the points and weights that are locally supported.
% Going descend direction (nu-uu) because that's how the book does
% And going ascend direction was giving some errors.
location = 0;
        for uu=0:pu
            location = location+1;
            pts(location,:) = P{nu-uu}(1:3);
            weight(location) = P{nu-uu}(4);
        end
% Fixing order to be ascending -> flipping up down


% Calculating the sum of weights
location = 0;
        for uu=0:pu
            location = location+1;
            % Total Weighting
            Q = Q + N(pu-uu+1)*weight(location);
            % Weight derivatives in regards to parametric coordinates
            dQ_du = dQ_du + dN(pu-uu+1)*weight(location);
        end

% Making the rational NURBS basis
location = 0;
        for uu=0:pu
            location = location+1;
            ratio = weight(location)/(Q*Q); % Rational NURBS Weight
            BASIS = N(pu-uu+1); % B-Spline Basis
            R(location) = BASIS*ratio*Q; % NURBS Basis
            % NURBS Basis derivatives with regards to parametric space
            dR_du(location,1) = (dN(pu-uu+1) - BASIS*dQ_du)*ratio;
         end

%% Jacobian
% From parameter space to physical space
location = 0;
for uu=0:pu
    location = location+1;
    dx_du = dx_du + P{nu-uu}(1)*dR_du(location,1);
end
du_dx = 1/dx_du; % Inverse of dx_du.

% Compute Derivatives of basis functions with respect to physical
% coordinates
dR(:,1) = dR_du*du_dx; % No for loop needed since this is 1D case.

Jacobian = dx_du; % Jacobian Matrix. Not really a matrix, since we're in 1D

J_mod = Jacobian;
J_2 = 0.5*(U(nu+1)-U(nu));
J = J_mod*J_2;
end
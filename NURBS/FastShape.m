%% INPUTS
% Knots - Cell Array containing the U, V and W Knot Vectors
% QPoints - Cell array containing the quadrature points
% QWeights - Cell array containing the quadrature rule weightings
%% OUTPUTS
% Basis - Cell array containing the NURBS univariate basis that support x
% Ders - Cell array containing the NURBS univ derivatives that support x
% DFinv - Cell array containing the DF^-1 inverse Jacobian matrix of each point
% detJ - Vector with the Jacobian Determinant multiplied by the QWeight.
function [Basis, Ders, DFinv, detJ] = FastShape(Model, Knots,QPoints,QWeights)
    Basis = cell(3,1);
    Ders = cell(3,1);
    DFinv = cell(3,1);
    detJ = cell(3,1);
    [q1, q2, q3] = size(QPoints);
    x = cell(q1,q2,q3);
    p = zeros(3,1);
    pu = Model.pu;
    pv = Model.pv;
    pw = Model.pw;
    p(1) = pu;
    p(2) = pv;
    p(3) = pw;
    P = Model.get_point_cell;
    for dim=1:3
        for ii=1:length(QPoints{dim})
            U = Knots{dim};
            qu = QPoints{dim}(ii);
            u = ((U(nu+1)-U(nu))*qu +U(nu+1)+U(nu))/2;
            pu = p(dim);
            n = length(U)-1-pu-1;
            span = FindSpanLinear(n, pu, u, U);
            nen = pu+1; % Number of Local Basis Functions
             % [-1,1] -> [0,1] remap
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
            J = J_mod*J_2*QWeights{dim}(ii);
            detJ{dim}(ii) = J;
            Basis{dim}(:,ii) = R;
            Ders{dim}(:,ii) = dR;
            DFinv{dim}(ii) = du_dx;
        end
    end
end
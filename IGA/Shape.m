function [ddR, Jacobian] = Shape(pu,u,U,P,MAX_ORDER_OF_DERIVATIVES)
%% INPUTS
% pu - The spline's degree
%  u - The point in [0,1] domain where the Shape funs are to be evaluated
%  U - Knot Vector in [0,1] range.
%  P - Control Points in U's univariate direction
% MAX_ORDER_OF_DERIVATIVES - The highest derivative to be computed.

%% OUTPUTS
% R - The shape functions
% dR - The first derivatives
% d2R - Second derivatives
% d3R - Third derivatives
% dnR - the n-th derivatives
% J - Jacobian

%% Algorithm
s = FindSpanLinear(length(U)-pu-2,pu,u,U); % Span in which lies u
Basis = DersBasisFun(s,u,pu,MAX_ORDER_OF_DERIVATIVES,U);

% Shape function, first derivatives and Jacobian

ACP = P(s-pu+1:s+1,:); % Active Control Points
Weights = ACP(:,4); % Weights of the Active Control Points
PP = ACP(:,1:3); 

clear ACP

W = Basis(1,:)*Weights; % Sum of Basis*Weights
R = Basis(1,:)/W; % Rational Basis Funs
x = Basis(1,:)*(PP.*Weights); % x y z coordinates

Jacobian = Basis(2,:)*(PP.*Weights);
J = norm(Jacobian);
dW = Basis(2,:)*Weights; % Sum of dN*Weights
dRdu = (1/(W^2))*(W*Basis(2,:) -dW*Basis(1,:));
dR = dRdu/J;
ddR = cell(MAX_ORDER_OF_DERIVATIVES,1);
ddR{1} = R;
ddR{2} = dR;
if MAX_ORDER_OF_DERIVATIVES > 1
    WW = [];
    WW(1) = W;
    WW(2) = dW;
    count = 1;
    while count < MAX_ORDER_OF_DERIVATIVES
        count = count +1;
        A = Weights'.*Basis(count+1,:);
        Wk = Basis(count+1,:)*Weights;
        WW(count+1) = Wk;
        sum = 0;
        for j=1:count
            sum = sum +nchoosek(count,j)*WW(j)*ddR{count-j+1};
        end
        ddR{count+1} = (A - sum)/(W*J^count);
    end
end
end
                






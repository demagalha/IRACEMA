function [ddR, J] = Shape(Model,u,MAX_ORDER_OF_DERIVATIVES)
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
switch Model.type
    case 'curve'
        p = Model.pu;
        U = {Model.U};
    case 'surf'
        p = [Model.pu; Model.pv];
        U = {Model.U; Model.V};
    case 'volume'
        p = [Model.pu; Model.pv; Model.pw];
        U = {Model.U; Model.V; Model.W};
end
for i=1:numel(p)
    s(i) = FindSpanLinear(length(U{i})-p(i)-2,p(i),u(i),U{i});
    Basis{i} = DersBasisFun(s(i),u(i),p(i),MAX_ORDER_OF_DERIVATIVES,U{i});
end
P = Model.get_point_cell;
switch Model.type
    case 'curve'
        ACP = P(s(1)-p(1)+1:s(1)+1);
    case 'surf'
        ACP = P(s(1)-p(1)+1:s(1)+1,s(2)-p(2)+1:s(2)+1);
    case 'volume'
        ACP = P(s(1)-p(1)+1:s(1)+1,s(2)-p(2)+1:s(2)+1,s(3)-p(3)+1:s(3)+1);
end
ACP = reshape(ACP,numel(ACP),1);
ACP = cell2mat(ACP); % Active Control Points
Weights = ACP(:,4); % Weights of the Active Control Points
PP = ACP(:,1:3); 

clear ACP
for i=1:MAX_ORDER_OF_DERIVATIVES+1
    access = length(u);
    b = Basis{access}(i,:);
    while access > 1
        b = kron(Basis{access-1}(i,:),b); % Kronecker Product Rules! :D
        access = access-1;
    end
    B(i,:) = b;
end
W = B(1,:)*Weights; % Sum of Basis*Weights
R = B(1,:)/W; % Rational Basis Funs
x = B(1,:)*(PP.*Weights); % x y z coordinates
Jacobian = B(2,:)*(PP.*Weights);
J = norm(Jacobian);
dW = B(2,:)*Weights; % Sum of dN*Weights
dRdu = (1/(W^2))*(W*B(2,:) -dW*B(1,:));
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
        A = Weights'.*B(count+1,:);
        Wk = B(count+1,:)*Weights;
        WW(count+1) = Wk;
        soma = 0;
        for j=1:count
            soma = soma +nchoosek(count,j)*WW(j)*ddR{count-j+1};
        end
        ddR{count+1} = (A - soma)/(W*J^count);
    end
end
end
                






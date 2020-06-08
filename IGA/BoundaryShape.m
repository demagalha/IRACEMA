function [R, dR, J] = BoundaryShape(Model,q,e,IEN,INN,BoundaryData)
% Shape function for Boundary Conditions
%% INPUTS
% Model - Geometry class object
% q - Gauss integration point in [-1,1]
% e - element number
% IEN - array that gives the global basis function number of a local basis
% from a determined element. GlobalBasisNumber = IEN(local,element)
% INN - array that gives the global basis function numbering given the
% unidimensional basis functions. It is essentially a sub2ind function.
% BoundaryData - array that carries Boundary information:
% BoundaryData(1) -> boundary parametric direction
% BoundaryData(2) -> integration parametric direction
% BoundaryData(3) -> boundary value. Either 0 or 1
%% OUTPUTS
% R - Basis functions
% dR - Derivative Basis functions
% J - Jacobian
%% Direction Arrays Declaration
% We start by making our lifes easier with BoundaryData information
bidx = BoundaryData(1);
pidx = BoundaryData(2);
bval = BoundaryData(3);
Knots = {Model.U, Model.V};
p = [Model.pu, Model.pv];
BKnot = Knots{bidx}; % Boundary Knot Vector
IKnot = Knots{pidx}; % Integration Knot Vector

%% Active Point Determination
% We start by determining where the support of the basis functions begin
n = INN(IEN(1,e),:);

ni = n(pidx);
pi = p(pidx);

nb = n(bidx);
pb = p(bidx);

BoundaryPointIndexes = nb+1:-1:nb-pb+1;
IntegrationPointIndexes = ni+1:-1:ni-pi+1;
if bidx == 1
    PointInterval = {BoundaryPointIndexes;IntegrationPointIndexes};
elseif bidx == 2
    PointInterval = {IntegrationPointIndexes;BoundaryPointIndexes};
else
    error('Invalid Boundary. This function is exclusive to 2D problems.');
end
P = Model.get_point_cell;
ActivePoints = cell((pi+1)*(pb+1),1);
for i = 1:numel(PointInterval{1})
    for j =1:numel(PointInterval{2})
        idx = sub2ind([numel(PointInterval{1}),numel(PointInterval{2})],i,j);
        ActivePoints{idx} = P{PointInterval{1}(i),PointInterval{2}(j)};
    end
end
ActivePoints = cell2mat(ActivePoints);



%% Shape Function Calculation
% Calculate the parametric direction integration point
ui = ((IKnot(ni+1) - IKnot(ni))*q +(IKnot(ni+1) + IKnot(ni)))*0.5;
ub = bval;

NI = DersBasisFun(ni-1,ui,pi,1,IKnot);
NB = zeros(pb+1,2);
if bval == 0
    NB(1,1) = 1;
    NB(2,1) = pb/(BKnot(pb+2)-BKnot(1));

end
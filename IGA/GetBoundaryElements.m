%% GetBoundaryElements
%  A function which returns the elements of the boundary of a 2D problem
% their parametric direction that is not being integrated, their side in
% the parametric direction and the lift value

%% INPUTS
% Model - Geometry class object
% lift_u - pair of values, corresponding to the lift at u = 0 and u = 1
% lift_v - pair of values, as above.

%% OUTPUTS
% BoundaryElements - n-by-4 matrix. Rows are elements of the boundary
% columns are:
% (:,1) - Element number in the IEN array
% (:,2) - Parametric Direction of the Boundary (1 for u, 2 for v)
% (:,3) - Value of the boundary. 0 for the u/v = 0 boundary, 1 for u/v = 1
% (:,4) - Lift value for the element, determined by lift_u and lift_v

function BoundaryElements = GetBoundaryElements(Model,direction,boundary,lift)
[INN, IEN, ~, ~] = Model.get_connectivity;

if boundary == 0
    u = find(INN(:,direction) == 1);
elseif boundary == 1
    u = find(INN(:,direction) == max(INN(:,direction)));
else
    error('boundary only accepts 0 or 1 as input value')
end
[~, idx, ~] = intersect(IEN,u);
[~, GammaElements] = ind2sub(size(IEN),idx);
GammaElements = unique(GammaElements);
boundary = repmat(boundary,length(GammaElements),1);
direction = repmat(direction,length(GammaElements),1);
lift = repmat(lift,length(GammaElements),1);
BoundaryElements = [GammaElements, direction, boundary,lift];

end
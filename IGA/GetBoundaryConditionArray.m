function BoundaryConditionArray = GetBoundaryConditionArray(Model,direction,boundary,lift)
[INN, IEN, ~, ~] = Model.get_connectivity;
ID = reshape(1:max(max(IEN)),1,max(max(IEN)));
if boundary == 0
    u = find(INN(:,direction) == 1);
elseif boundary == 1
    u = find(INN(:,direction) == max(INN(:,direction)));
else
    error('boundary only accepts 0 or 1 as input value')
end
lift = repmat(lift,length(u),1);
BoundaryConditionArray = [u, lift];
end
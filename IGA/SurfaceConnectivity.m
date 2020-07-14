function [global_basis_index, element_local_mapping, element_range] = ...
    SurfaceConnectivity(GeometryObj)
array_size = size(GeometryObj.PX);
ndof = prod(array_size);
[n(:), m(:)] = ind2sub(array_size,1:ndof);
global_basis_index = [n',m']; % global basis index == INC vector

p = GeometryObj.PolynomialOrder;
ELEMENT_DOFS = prod(p+1);
Knots =  GeometryObj.KnotVectorCell;
elements_per_direction = zeros(size(p));
basis_spans = cell(size(Knots));
element_ranges = basis_spans;
for i=1:length(Knots)
    elements_per_direction(i) = length(unique(Knots{i}))-1;
    [element_ranges{i}, basis_spans{i}] = KnotConnectivity(p(i),Knots{i}); 
end
u_spans = cell2mat(basis_spans(1));
u_ranges = cell2mat(element_ranges(1));
v_spans = cell2mat(basis_spans(2));
v_ranges = cell2mat(element_ranges(2));
ELEMENTS = prod(elements_per_direction);
element_local_mapping = zeros(ELEMENT_DOFS, ELEMENTS);
element_range = zeros(ELEMENTS,2,2);
e_out = cell(size(Knots));
[e_out{:}] = ind2sub(elements_per_direction,1:ELEMENTS);
e_out = cell2mat(e_out');
e_out = e_out'; % Global element number from parametric element nml
for e=1:ELEMENTS
    tmp = e_out(e,:);
    ei = tmp(1);
    ej = tmp(2);
    element_u_spans = u_spans(ei,:)';
    element_u_ranges = u_ranges(ei,:);
    element_v_spans = v_spans(ej,:)';
    element_v_ranges = v_ranges(ej,:);
    column_1 = repmat(element_u_spans,length(element_v_spans),1);
    column_2 = repmat(element_v_spans,length(element_u_spans),1);
    column_2 = sort(column_2,'asc');
    column_3 = [column_1, column_2];
    element_spans = column_3;
    element_basis = sub2ind(array_size,element_spans(:,1), ... 
        element_spans(:,2));
    element_local_mapping(:,e) = element_basis;
    element_range(e,:,1) = element_u_ranges;
    element_range(e,:,2) = element_v_ranges;
end
end
function [element_ranges, element_connectivity] = KnotConnectivity(P, ...
    knot_vector)
NUMBER_OF_ELEMENTS = length(unique(knot_vector))-1;
element_ranges = zeros(NUMBER_OF_ELEMENTS,2);
element_knot_indices = zeros(NUMBER_OF_ELEMENTS,2);
element_connectivity = zeros(NUMBER_OF_ELEMENTS,P+1);
element = 1;
PREVIOUS_VALUE = 0;
for i=1:length(knot_vector)
    CURRENT_VALUE = knot_vector(i);
    if knot_vector(i) ~= PREVIOUS_VALUE
        element_ranges(element,:) = [PREVIOUS_VALUE CURRENT_VALUE];
        element_knot_indices(element,:) = [i-1 i];
        element = element +1;
    end
    PREVIOUS_VALUE = CURRENT_VALUE;
end

REPEATED_KNOTS = 0;
for e=1:NUMBER_OF_ELEMENTS
    indices = element_knot_indices(e,1)-P+1:element_knot_indices(e,1);
    previous_values = knot_vector(indices);
    current_values = ones(1,P)*knot_vector(element_knot_indices(e,1));
    if isequal(previous_values,current_values) && ...
            length(nonzeros(previous_values))>1
        REPEATED_KNOTS = REPEATED_KNOTS +1;
    end
    element_connectivity(e,:) = ... 
        element_knot_indices(e,1)-P:element_knot_indices(e,1);
end
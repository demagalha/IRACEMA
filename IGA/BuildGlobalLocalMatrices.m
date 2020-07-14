function [id, lm] = ...
    BuildGlobalLocalMatrices(element_local_mapping,NUMBER_OF_SPACE_D)
    [nen, nel] = size(element_local_mapping);
    MAX_BASIS = max(max(element_local_mapping));
    id = reshape(1:MAX_BASIS*NUMBER_OF_SPACE_D,NUMBER_OF_SPACE_D,MAX_BASIS);
    lm = zeros(NUMBER_OF_SPACE_D*nen,nel);
    for i=1:nel
        lm(:,i) = reshape(id(:,element_local_mapping(:,i)), ...
            NUMBER_OF_SPACE_D*nen,1);
    end
end
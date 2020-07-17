function [d, F] = ApplyDirichletBCs(K,F,BoundaryBasis,boundary_values)
d = zeros(size(F));
FreeDOFS = setdiff(1:length(d),BoundaryBasis);
d(BoundaryBasis) = boundary_values;
F(FreeDOFS) = F(FreeDOFS) - K(FreeDOFS,BoundaryBasis)*boundary_values;
d(FreeDOFS) = K(FreeDOFS,FreeDOFS)\F(FreeDOFS);
end
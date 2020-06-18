function [d, F] = DirichletBC(K,F,DirichletBoundaryArray)
d = zeros(size(F));

BoundaryDOFS = DirichletBoundaryArray(:,1);
FreeDOFS = setdiff(1:length(d),BoundaryDOFS);

d(BoundaryDOFS) = DirichletBoundaryArray(:,2);
F(FreeDOFS) = F(FreeDOFS) -K(FreeDOFS,BoundaryDOFS)*DirichletBoundaryArray(:,2);
d(FreeDOFS) = K(FreeDOFS,FreeDOFS)\F(FreeDOFS);
end
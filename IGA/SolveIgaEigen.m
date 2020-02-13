function [autovector,omega] = SolveIgaEigen(K,M,bc,nmode)

%K global stiffness matrix, ndof x ndof
%M global mass matrix, ndof x ndof
%bc, vector of bcs x 2, first bc(i,1) is the index of de dof and bc(i,2) =
%is the prescribed displacement %%currently just bcs, without prescribed
%values
%nmode number of modes

N_DOF = size(K,1);

bc = sort(bc,'descend');


%remove columns and lines related to the dof
for i=1:numel(bc)
    K(:,bc(i)) = [];
    K(bc(i),:) = [];
    M(:,bc(i)) = [];
    M(bc(i),:) = [];
end

[autovector,omega] = eigs(K,M,nmode,'sm');
omega = sqrt(diag(omega));
bc = sort(bc,'ascend');



autovector = BoundariesPostProcess(autovector,bc);

end

    

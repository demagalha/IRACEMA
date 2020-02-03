clear all
close all
clc

%% Pre-processing
    load timoshenkoruled_coarse.mat
    % hpk refine the model ------------------------------------------------
    Model.DegreeElevate(4,1);
    Model.KnotRefine([0.1:0.1:0.9],1);
    % Material Properties -------------------------------------------------
    YOUNG_MODULUS = 210*10^9;
    RHO = 7860; 
    POISSON = 0.3;
    D = get_matprop_matrix(1,YOUNG_MODULUS,POISSON); % Isotropic D Matrix
%% Assembly
    [K, M, IEN] = Assemble(Model,D,RHO);

%% Boundary Conditions
    bc_clamp_y = 0; % Beam is clamped in y=0
%     bc_sup_y1 = 1; % Beam is simply supported in y=1 and z=0
%     bc_sup_y2 = 0; % Beam is simply supported in y=0 and z=0
    bc_sup_z = 0; % Beam is simply supported in z=0 and y=1
    P = Model.get_point_cell;
    clamped_pts = [];
    for i=1:numel(P)
        if (P{i}(2) == bc_clamp_y)
            clamped_pts = [clamped_pts i];
        end
    end
    ID = reshape(1:max(max(IEN))*3,3,max(max(IEN)));
    clamp_ctrl_pts = reshape(ID(:,clamped_pts),numel(ID(:,clamped_pts)),1);
%     support_ctrl_pts = reshape(ID(3,supported_pts),numel(ID(3,supported_pts)),1);
    bc = [clamp_ctrl_pts]; % boundary conditions
    bc = sort(bc,'descend');
    K = full(K); % Easier to manipulate in MATLAB than sparse
    for i=1:numel(bc)
        K(bc(i),bc(i)) = 1e30;
    end
    K = sparse(K);
    % Quick comment on BC: I don't like manipulating big matrice sizes and
    % etc, as it costs too much computational power. With this way of
    % enforcing BCs, since we're only dealing with fixing/clamping style,
    % adding a very high rigidity to the DOF in question is faster/more
    % elegant.
%% Post-Processing
    [autovector,ome] = eigs(K,M,10,'sm');
    omega = sqrt(ome);
    freq = omega/(2*pi);
    clearvars -except K M autovector freq omega Model ID
    save('Example_Beam3D.mat')
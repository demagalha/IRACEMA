clear all
close all
clc

%% Pre-processing

    walue = zeros(7,20);
    for r=1:20
    tic;
    load RectangularLeissa15.mat
    % hpk refine the model ------------------------------------------------
    Model.DegreeElevate(4,1);
    Model.DegreeElevate(4,2);
    Model.KnotRefine(1/(r+1):1/(r+1):1-1/(r+1),1);
    Model.KnotRefine(1/(r+1):1/(r+1):1-1/(r+1),2);
    % Material Properties -------------------------------------------------
    YOUNG_MODULUS = 210*10^9;
    RHO = 7860; 
    POISSON = 0.3;
    D = get_matprop_matrix(1,YOUNG_MODULUS,POISSON); % Isotropic D Matrix
%% Assembly
    [K, M, IEN] = Assemble(Model,D,RHO);

%% Boundary Conditions
    clamp_x = [0 1];
    clamp_y = [0 1.5];
    clamp = [clamp_x;clamp_y];
    P = Model.get_point_cell;
    clamped_pts = [];
    for i=1:numel(P)
        if (abs(P{i}(1) - clamp_x(1)) < sqrt(eps)) || (abs(P{i}(1) - clamp_x(2)) < sqrt(eps))
            clamped_pts = [clamped_pts i];
        elseif  (abs(P{i}(2) - clamp_y(1)) < sqrt(eps)) || (abs(P{i}(2) - clamp_y(2)) < sqrt(eps))
            clamped_pts = [clamped_pts i];
        end
    end
    ID = reshape(1:max(max(IEN))*3,3,max(max(IEN)));
    clamp_ctrl_pts = reshape(ID(:,clamped_pts),numel(ID(:,clamped_pts)),1);
%     support_ctrl_pts = reshape(ID(3,supported_pts),numel(ID(3,supported_pts)),1);
    bc = [clamp_ctrl_pts]; % boundary conditions
    K = full(K); % Easier to manipulate in MATLAB than sparse
    M = full(M);
    for i=1:numel(bc)
        K(bc(i),bc(i)) = 1e30;
    end
    K = sparse(K);
    M = sparse(M);
    % Quick comment on BC: I don't like manipulating big matrice sizes and
    % etc, as it costs too much computational power. With this way of
    % enforcing BCs, since we're only dealing with fixing/clamping style,
    % adding a very high rigidity to the DOF in question is faster/more
    % elegant.
%% Post-Processing
    [autovector,ome] = eigs(K,M,10,'sm');
    bc = sort(bc,'ascend');
%     for i=1:numel(bc)
%         autovector = [auto
%     end
    omega = diag(sqrt(ome));
    freq = omega/(2*pi);
    h = 0.02;
    rho = RHO*h;
    D = YOUNG_MODULUS*h*h*h/(12*(1-POISSON*POISSON));
    walue([1:6],r) = omega([1:6])*sqrt(rho/D);
    [~, ~, nel, ~] = Model.get_connectivity;
    walue(6,r) = toc;
    walue(7,r) = nel;
    walue(:,r)
    end
    clearvars -except K M autovector freq omega Model ID walue nel
    save('Example_SquarePlateP5.mat')
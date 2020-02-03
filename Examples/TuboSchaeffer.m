clear all
close all
clc

%% Preprocessing
load cano_schaeffer.mat

% Refinar o modelo
Model.DegreeElevate(1,1)
Model.DegreeElevate(2,3);
Model.DegreeElevate(2,2);
Model.KnotRefine([0.01:0.01:0.99],3);
[INN, ~, ~, ~] = Model.get_connectivity;
% Material Properties
    YOUNG_MODULUS = 2e11;
    RHO = 7850; 
    POISSON = 0.3;
    D = get_matprop_matrix(1,YOUNG_MODULUS,POISSON); % Isotropic Deformation
%% Assembly
        [K, M, IEN] = Assemble(Model,D,RHO);
%% Apply Forces
    P = Model.get_point_cell;
    F = zeros(numel(INN),1);
    ID = reshape(1:max(max(IEN))*3,3,max(max(IEN)));
    constNod = [];
    for i=1:103
        constNod = [constNod; [40:52]'+(i-1)*52];
    end
    constNod = reshape(constNod,numel(constNod),1);
    ForceValues = zeros(numel(constNod),2);
    for i=1:numel(constNod)
        angle = atan(P{constNod(i)}(2)/P{constNod(i)}(1));
        ForceValues(i,:) = [cos(angle), sin(angle)];
    end
    ForcesLocationX = reshape(ID(1,constNod),numel(ID(1,constNod)),1);
    ForcesLocationY = reshape(ID(2,constNod),numel(ID(2,constNod)),1);
    for  i=1:numel(constNod)
        F(ForcesLocationX(i)) = ForceValues(i,1);
        F(ForcesLocationY(i)) = ForceValues(i,2);
    end
        

%% Post-Processing
    K = sparse(K);
    M = sparse(M);
    F = sparse(F);
    freq = 30:10:5000;
    x = zeros(numel(INN),length(freq));
    fprintf('Começando Inversão Matricial \n')
    for i=1:length(freq)
        omega = 2*pi*freq(i);
        ComplexStiffness = K -(omega^2)*M;
        x(:,i) = F\ComplexStiffness;
    end
    fprintf('Inversão Realizada \n')
  
%     Find Deflection Places
    PointsOfInterest = [];
    for i=1:103
        PointsOfInterest = [PointsOfInterest; [4;10]+(i-1)*52];
    end
    
    clearvars -except K M x freq Model ID PointsOfInterest
    save('Ricardo.mat')
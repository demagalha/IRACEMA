function Model = BuildHoledPlate
% P{1} = [-1 0 0 1];
% P{2} = [-2.5 0 0 1];
% P{3} = [0 L 0 1];
% C1 = Geometry('curve',1,[0 0 0.5 1 1],P);
% 
% B{1} = [-R 0 0 1];
% B{2} = [-R R 0 sqrt(2)/2];
% B{3} = [0 R 0 1];
% C2 = Geometry('curve',2,[0 0 0 1 1 1],B);
% 
% Model = geo_ruled(C1,C2,C1,C1);
controlPts        = zeros(4,3,2);

controlPts(1,1,:) = [-1 0];
controlPts(1,2,:) = [-2.5 0];
controlPts(1,3,:) = [-4 0];

controlPts(2,1,:) = [-1 0.4142135623730951];
controlPts(2,2,:) = [-2.5 0.75];
controlPts(2,3,:) = [-4 4];

controlPts(3,1,:) = [-0.4142135623730951 1];
controlPts(3,2,:) = [-0.75 2.5];
controlPts(3,3,:) = [-4 4];

controlPts(4,1,:) = [0 1];
controlPts(4,2,:) = [0 2.5];
controlPts(4,3,:) = [0 4];

uKnot = [0 0 0 0.5 1 1 1];
vKnot = [0 0 0 1 1 1];

noPtsX = 4;
noPtsY = 3;

controlPts = reshape(controlPts,12,2);

p     = 2;
q     = 2;

cont = 0.5*(1+1/sqrt(2));

weights = [1; cont; cont; 1;
           1; 1; 1; 1;
           1; 1; 1; 1];
P = [controlPts, zeros(length(controlPts),1), weights];
P = num2cell(P,2);
P = reshape(P,noPtsX,noPtsY);
Model = Geometry('surf',p,uKnot,q,vKnot,P);
end
function [R, d2R, J] = BeamShape(Model,qu,element,P,IEN,INN)

%% Setup
pu = Model.pu;
nen = (pu+1);
tmp = INN(IEN(1,element),:);
nu = tmp(1);
clear tmp
U = Model.U;
u = ((U(nu+1)-U(nu))*qu +U(nu+1)+U(nu))/2; % Change of variables
R = zeros(nen,1);
d2R = zeros(nen,1);
Q = 0;
dQdu = 0;
d2Qdu2 = 0;
weight = zeros(1,nen);
dxdu = 0;

%% B-Spline Basis Functions
NU = DersBasisFun(nu-1,u,pu,2,U);
N = NU(1,:);
dN = NU(2,:);
d2N = NU(3,:);
clear NU

%% NURBS Basis Functions

location = 0;
for uu=0:pu
    location = location+1;
    weight(location) = P{nu-uu}(4);
end

location = 0;
for uu=0:pu
    location = location+1;
    Q = Q+ N(pu-uu+1)*weight(location);
    dQdu = dQdu + dN(pu-uu+1)*weight(location);
    d2Qdu2 = d2Qdu2 +d2N(pu-uu+1)*weight(location);
end

location = 0;
for uu=0:pu
    location = location+1;
    ww = weight(location)/(Q*Q);
    R(location) = N(pu-uu+1)*ww*Q; % NURBS Basis
    dR(location,1) = (dN(pu-uu+1) - N(pu-uu+1)*dQdu)*ww;
    d2R(location,1) = (1/Q)*(weight(location)*d2N(pu-uu+1) -2*dR(location)*dQdu -R(location)*d2Qdu2);
end

%% Jacobian
location = 0;
for uu=0:pu
    location = location+1;
    dxdu = dxdu +P{nu-uu}(1)*dR(location);
end
dudx = inv(dxdu);

d2R(:,1) = d2R*(dudx)^2;
Jmod = dxdu;
J2 = 0.5*(U(nu+1)-U(nu));
J = Jmod*J2;
end

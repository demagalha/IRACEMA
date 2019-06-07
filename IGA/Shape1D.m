function [R, dR_dx, J] = Shape1D(xi_tilde,e,p,P,KV_Xi,INN,IEN)
 
nen = (p+1) ; % Number of Local Shape Functions
%
R(nen) = 0; %array of trivariate NURBS basis funs
dR_dx(nen) = 0; %derivative
 
J = 0; %jacobian determinant
%
 
ni = 0;
xi = 0; % parametric coordinates
N(p+1) = 0; %array of univariate B-spline basis funs
 
dN_dxi(p+1) = 0; %univariate B-spline
    %function derivatives w.r.t. %appropriate parametric coordinates
dR_dxi(nen) = 0; %trivariate NURBS function derivatives
                    %w.r.t. parametric coordinates
dx_dxi = 0; %derivative of physical coordinates
                %w.r.t. parametric coordinates
               
dxi_dx = 0; %inverse of dx_dxi
dxi_dtildexi = 0; %derivative of parametric coordinates
                       %w.r.t parent element coordinates
J_mat = 0; %jacobian "matrix"
i = 0;
loc_num = 0 ; %local basis function counter
sum_xi = 0; sum_eta = 0; sum_zeta = 0; sum_tot = 0; %dummy sums
 
 
%NURBS coordinates; see algo 7(book)
%ni = INN(IEN(e,1),1);
%nj = INN(IEN(e,1),2);
%nk = INN(IEN(e,1),3);
 
ni = INN(IEN(1,e),1);

 
%calculate parametric coordinates from parent element coordinates
%knot vectors KV_Xi, KV_Eta and KV_Zeta and
%parent element coordinates xi_tilde, eta_tilde, zeta_tilde
%are given as input
 
xi = ((KV_Xi(ni+1)-KV_Xi(ni))*xi_tilde + (KV_Xi(ni+1) + KV_Xi(ni)))/2;

 
 
%%%part 2
 
 
N1 = DersBasisFun(ni-1,xi,p,1,KV_Xi);
N = N1(1,:);
dN_dxi = N1(2,:);
 
clear N1
 
        for i=0:p
            loc_num = loc_num +1;
            R(loc_num) = N(p+1-i) * P{ni-i,1,1}(4);
            sum_tot = sum_tot + R(loc_num);
            dR_dxi(loc_num,1) = dN_dxi(p+1)*P{ni-i}(4);
            sum_xi = sum_xi + dR_dxi(loc_num,1);
        end 
%divide by denomominators to complete definitions of function
% and derivatives w.r.t. parametric coordinates
 
for loc_num = 1 : nen
    R(loc_num) = R(loc_num)/sum_tot;
    dR_dxi(loc_num) = (dR_dxi(loc_num)*sum_tot - R(loc_num)*sum_xi)/ (sum_tot * sum_tot);
end
 
 
%%%part 3
 
loc_num = 0;
 
        for i=0:p
            loc_num = loc_num+1;
            dx_dxi = dx_dxi + P{ni-i}(1)*dR_dxi(loc_num);
        end
 
%dx_dxi
dxi_dx = 1/(dx_dxi);
 
%compute derivatives of basis functions
%with respect to physical coordinates
 
for loc_num=1:nen
    dR_dx(loc_num) = dR_dx(loc_num) + dR_dxi(loc_num)*dxi_dx;
end
 
 
%gradient of mapping from parent element to parameter space
dxi_dtildexi = (KV_Xi(ni+1)-KV_Xi(ni))/2;

J_mat = J_mat + dx_dxi*dxi_dtildexi;
J = abs(J_mat);
end
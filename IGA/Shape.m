function [R, dR_dx, J] = Shape(xi_tilde, eta_tilde, zeta_tilde,e,p,q,r,P,KV_Xi,KV_Eta,KV_Zeta,INN,IEN)
 
nen = (p+1)*(q+1)*(r+1); % Number of Local Shape Functions
%
R(nen,1) = 0; %array of trivariate NURBS basis funs
dR_dx(nen,3) = 0; %derivative
 
J = 0; %jacobian determinant
%
 
ni = 0; nj = 0; nk = 0;
xi = 0; eta = 0; zeta = 0; % parametric coordinates
N(p+1) = 0; M(q+1) = 0; L(r+1) = 0; %array of univariate B-spline basis funs
 
dN_dxi(p+1) = 0; %univariate B-spline
dM_deta(q+1) = 0; %function derivatives w.r.t.
dL_dzeta(r+1) = 0; %appropriate parametric coordinates
dR_dxi(nen,3) = 0; %trivariate NURBS function derivatives
                    %w.r.t. parametric coordinates
dx_dxi(3,3) = 0; %derivative ofp hysical coordinates
                %w.r.t. parametric coordinates
               
dxi_dx(3,3) = 0; %inverse of dx_dxi
dxi_dtildexi(3,3) = 0; %derivative of parametric coordinates
                       %w.r.t parent element coordinates
J_mat(3,3) = 0; %jacobian matrix
i = 0; j = 0; k = 0; aa = 0; bb = 0; cc = 0;
loc_num = 0 ; %local basis function counter
sum_xi = 0; sum_eta = 0; sum_zeta = 0; sum_tot = 0; %dummy sums
 
 
%NURBS coordinates; see algo 7(book)
%ni = INN(IEN(e,1),1);
%nj = INN(IEN(e,1),2);
%nk = INN(IEN(e,1),3);
 
ni = INN(IEN(1,e),1);
nj = INN(IEN(1,e),2);
nk = INN(IEN(1,e),3);
 
 
%calculate parametric coordinates from parent element coordinates
%knot vectors KV_Xi, KV_Eta and KV_Zeta and
%parent element coordinates xi_tilde, eta_tilde, zeta_tilde
%are given as input
 
xi = ((KV_Xi(ni+1)-KV_Xi(ni))*xi_tilde + (KV_Xi(ni+1) + KV_Xi(ni)))/2;
eta = ((KV_Eta(nj+1)-KV_Eta(nj))*eta_tilde + (KV_Eta(nj+1) + KV_Eta(nj)))/2;
zeta = ((KV_Zeta(nk+1)-KV_Zeta(nk))*zeta_tilde + (KV_Zeta(nk+1) + KV_Zeta(nk)))/2;
 
 
 
 
%%%part 2
 
 
N1 = DersBasisFun(ni-1,xi,p,1,KV_Xi);
M1 = DersBasisFun(nj-1,eta,q,1,KV_Eta);
L1 = DersBasisFun(nk-1,zeta,r,1,KV_Zeta);
N = N1(1,:);
dN_dxi = N1(2,:);
M = M1(1,:);
dM_deta = M1(2,:);
L = L1(1,:);
dL_dzeta = L1(2,:);
 
clear N1 M1 L1
 
for k=0:r
    for j=0:q
        for i=0:p
            loc_num = loc_num +1;
           
            %fprintf('p+1-i : %d , q+1-j : %d , r+1-k %d \t ni-i : %d , nj-j : %d , nk-k : %d\n', p+1-i, q+1-j, r+1-k, ni-i, nj-j,nk-k);
            R(loc_num) = N(p+1-i)*M(q+1-j)*L(r+1-k) * P{ni-i,nj-j,nk-k}(4);
            sum_tot = sum_tot + R(loc_num);
           
            dR_dxi(loc_num,1) = dN_dxi(p+1-i)*M(q+1-j)*L(r+1-k) * P{ni-i,nj-j,nk-k}(4);
            sum_xi = sum_xi + dR_dxi(loc_num,1);
           
            dR_dxi(loc_num,2) = N(p+1-i)*dM_deta(q+1-j)*L(r+1-k) * P{ni-i,nj-j,nk-k}(4);
            sum_eta = sum_eta + dR_dxi(loc_num,2);
           
            dR_dxi(loc_num,3) = N(p+1-i)*M(q+1-j)*dL_dzeta(r+1-k) * P{ni-i,nj-j,nk-k}(4);
            sum_zeta = sum_zeta + dR_dxi(loc_num,3);
        end
    end
end
 
%divide by denomominators to complete definitions of function
% and derivatives w.r.t. parametric coordinates
 
for loc_num = 1 : nen
    R(loc_num) = R(loc_num)/sum_tot;
   
    dR_dxi(loc_num,1) = (dR_dxi(loc_num,1)*sum_tot - R(loc_num)*sum_xi)/ (sum_tot * sum_tot);
    dR_dxi(loc_num,2) = (dR_dxi(loc_num,2)*sum_tot - R(loc_num)*sum_eta)/ (sum_tot * sum_tot);
    dR_dxi(loc_num,3) = (dR_dxi(loc_num,3)*sum_tot - R(loc_num)*sum_zeta)/ (sum_tot * sum_tot);
end
 
 
%%%part 3
 
loc_num = 0;
 
for k=0:r
    for j=0:q
        for i=0:p
            loc_num = loc_num+1;
           
            for aa=1:3
                for bb=1:3
                    dx_dxi(aa,bb) = dx_dxi(aa,bb) + P{ni-i,nj-j,nk-k}(aa) * dR_dxi(loc_num,bb);
                end
            end
        end
    end
end
 
%dx_dxi
dxi_dx = inv(dx_dxi);
 
%compute derivatives of basis functions
%with respect to physical coordinates
 
for loc_num=1:nen
    for aa=1:3
        for bb=1:3
            dR_dx(loc_num,aa) = dR_dx(loc_num,aa) + dR_dxi(loc_num,bb)*dxi_dx(bb,aa);
        end
    end
end
 
 
%gradient of mapping from parent element to parameter space
dxi_dtildexi(1,1) = (KV_Xi(ni+1)-KV_Xi(ni))/2;
dxi_dtildexi(2,2) = (KV_Eta(nj+1)-KV_Eta(nj))/2;
dxi_dtildexi(3,3) = (KV_Zeta(nk+1)-KV_Zeta(nk))/2;
 
for aa=1:3
    for bb=1:3
        for cc=1:3
            J_mat(aa,bb) = J_mat(aa,bb) + dx_dxi(aa,cc)*dxi_dtildexi(cc,bb);
        end
    end
end
 
 
J = det(J_mat);
 
end
                






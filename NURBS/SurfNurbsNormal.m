function N = SurfNurbsNormal(Model,u,v)

%given a parameter (u,v) returns the normal to the surface at the point specified (u,v)

wx = zeros(Model.nu + 1,Model.nv + 1);
wy = wx;
wz = wx;
w = wx;
for i=1:Model.nu + 1
	for j=1:Model.nv + 1
		wx(i,j) = Model.Pw(i,j).x;
		wy(i,j) = Model.Pw(i,j).y;
		wz(i,j) = Model.Pw(i,j).z;
		w(i,j) = Model.Pw(i,j).w;
	end
end

Aders.x = SurfDerivAlg1(Model.nu,Model.pu,Model.U,Model.nv,Model.pv,Model.V,wx,u,v,1);
Aders.y = SurfDerivAlg1(Model.nu,Model.pu,Model.U,Model.nv,Model.pv,Model.V,wy,u,v,1);
Aders.z = SurfDerivAlg1(Model.nu,Model.pu,Model.U,Model.nv,Model.pv,Model.V,wz,u,v,1);
wders = SurfDerivAlg1(Model.nu,Model.pu,Model.U,Model.nv,Model.pv,Model.V,w,u,v,1);

D.x = RatSurfDerivs(Aders.x,wders,1);
D.y = RatSurfDerivs(Aders.y,wders,1);
D.z = RatSurfDerivs(Aders.z,wders,1);

Su = [D.x(2,1) D.y(2,1) D.z(2,1)];
Sv = [D.x(1,2) D.y(1,2) D.z(1,2)];

N = cross(Su,Sv);
N = N/norm(N);
end 
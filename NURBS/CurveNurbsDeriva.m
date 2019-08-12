function SS = CurveNurbsDeriva(Model,u,d)

%given a parameter (u) returns the d derivative [x,y,z]


wx = zeros(1,numel(Model.Pw));
wy = wx;
wz = wx;
w = wx;
for i=1:numel(Model.Pw)
	wx(i) = Model.Pw(i).x;
	wy(i) = Model.Pw(i).y;
	wz(i) = Model.Pw(i).z;
	w(i) = Model.Pw(i).w;
end 

Aders.x = CurveDerivsAlg1dot2(Model.nu,Model.pu,Model.U,wx,u,d);
Aders.y = CurveDerivsAlg1dot2(Model.nu,Model.pu,Model.U,wy,u,d);
Aders.z = CurveDerivsAlg1dot2(Model.nu,Model.pu,Model.U,wz,u,d);
wders = CurveDerivsAlg1dot2(Model.nu,Model.pu,Model.U,w,u,d);



S.x = RatCurveDerivs2(Aders.x,wders,d);
S.y = RatCurveDerivs2(Aders.y,wders,d);
S.z = RatCurveDerivs2(Aders.z,wders,d);

SS = [S.x(d+1) S.y(d+1) S.z(d+1)]; 

end

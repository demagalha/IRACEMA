function normal = SurfaceNormal(Model,point)

  S = SurfDerivAlg1(Model.n,Model.pu,Model.U,Model.nv,Model.pv,Model.V,PX,point(1),point(2),1);
  S2 = SurfDerivAlg1(Model.n,Model.pu,Model.U,Model.nv,Model.pv,Model.V,Model.PY,point(1),point(2),1);
  S3 = SurfDerivAlg1(Model.n,Model.pu,Model.U,Model.nv,Model.pv,Model.V,Model.PZ,point(1),point(2),1);
  
 U = [S(2,1) S2(2,1) S3(2,1)];
 V = [S(1,2) S2(1,2) S3(1,2)];
 
 normal = cross(U,V);
end
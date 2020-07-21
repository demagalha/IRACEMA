function N = CurveNurbsNormal(Model,u)

switch Model.type

	case 'curve'
		r1 = CurveNurbsDeriva(Model,u,1);
		r2 = CurveNurbsDeriva(Model,u,2);
		T1 = cross(r1,cross(r2,r1));
		N = T1/norm(T1);
	otherwise
		disp('Input model must be a curve');
end

end

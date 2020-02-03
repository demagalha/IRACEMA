function revol = geo_revolve(Model)

theta = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4, 2*pi];

scaling = [0 1 0 1 0 1 0 1 0];
scaling_factor = {[1 1 1], [sqrt(2) 1 1], [1 1 1], [sqrt(2) 1 1], [1 1 1], [sqrt(2) 1 1], [1 1 1], [sqrt(2) 1 1] [1 1 1]};

	if strcmp(Model.type,'surf')
		S = Model.get_point_cell;
        pu = Model.pu;
		U = Model.U;
        pv = Model.pv;
		V = Model.V;
	
        ref = S{2,2};
        
		circle = geo_circle([0 0 0],0,'xy');
		C = circle.get_point_cell;
        W = circle.U;
		
		tam = size(S);
		P = cell(tam(1),tam(2),9);
     
        
		for i=1:tam(1)
			for j=1:tam(2)
				for k=1:9
                    
                    if scaling(k) == 1
                       Model2 = geo_scaling(Model,scaling_factor{k});
                    else
                        Model2 = Model;
                    end
                    rot = geo_rotate(Model2,[0 0 1], theta(k));
                    SS = rot.get_point_cell;
					P{i,j,k} = [SS{i,j}(1) + C{k}(1) ; SS{i,j}(2) + C{k}(2) ; SS{i,j}(3) + C{k}(3) ; SS{i,j}(4) * C{k}(4)]; 
                end
            end
        end
        
        revol = Geometry('volume',pu,U,pv,V,2,W,P);
    end
end

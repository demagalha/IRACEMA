classdef Geometry < handle
    properties
        type;
        nu;
        pu;
        U;
        nv;
        pv;
        V;
        nw;
        pw;
        W;
        PX;
        PY;
        PZ;
        weight;
        Pw;
        PlotProp;
    end
    
    methods
        
        function obj = Geometry(varargin)
            
            if nargin == 0
                return
            end
            
            switch varargin{1}
				case 'curve'
					obj.type = varargin{1};
					obj.pu = varargin{2};
					obj.U =  varargin{3};

					[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdatePTS(obj, varargin{4});
					obj.nu = size(varargin{4},2)-1;

					temp(1:obj.nu +1) = CPOINT(0,0,0,0,1);
					for i=1:obj.nu+1
						temp(i) = CPOINT(obj.PX(i),obj.PY(i),obj.PZ(i),obj.weight(i),0);
					end
					obj.Pw = temp;
					obj.NormalizeKnotVector;
					obj.PlotProp = struct('RGB',[0 0 1],'LineSize',2,'ControlRGB',[0 0 0],'IsoRGB',[],'MarkerCPT','o','MarkerCPTRGB',[1 0 0]);
					
				case 'surf'
					obj.type = varargin{1};
					obj.pu = varargin{2};
					obj.U = varargin{3};
					obj.pv = varargin{4};
					obj.V = varargin{5};

					[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdatePTS(obj, varargin{6});
					obj.nu = size(varargin{6},1)-1;
					obj.nv = size(varargin{6},2)-1;

					temp(1:obj.nu+1,1:obj.nv+1) = CPOINT(0,0,0,0,1);
					for i=1:obj.nu+1
						for j=1:obj.nv+1
							temp(i,j) = CPOINT(obj.PX(i,j),obj.PY(i,j),obj.PZ(i,j),obj.weight(i,j),0);
						end
					end
					obj.Pw = temp;
					obj.NormalizeKnotVector;
					obj.PlotProp = struct('RGB',[.1 .9 .1],'LineSize',2,'ControlRGB',[1 0 0],'IsoRGB',[0 0 0],'MarkerCPT','o','MarkerCPTRGB',[1 0 0]);
					
				case 'volume'
					obj.type = varargin{1};
					obj.pu = varargin{2};
					obj.U = varargin{3};
					obj.pv = varargin{4};
					obj.V = varargin{5};
					obj.pw = varargin{6};
					obj.W = varargin{7};

					[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdatePTS(obj, varargin{8});
					obj.nu = size(varargin{8},1)-1;
					obj.nv = size(varargin{8},2)-1;
					obj.nw = size(varargin{8},3)-1;

					temp(1:obj.nu+1,1:obj.nv+1,1:obj.nw+1) = CPOINT(0,0,0,0,1);
					for i=1:obj.nu+1
						for j=1:obj.nv+1
							for k=1:obj.nw+1
								temp(i,j,k) = CPOINT(obj.PX(i,j,k),obj.PY(i,j,k),obj.PZ(i,j,k),obj.weight(i,j,k),0);
							end
						end
					end
					obj.Pw = temp;
					obj.NormalizeKnotVector;
					obj.PlotProp = struct('RGB',[.1 .9 .1],'LineSize',2,'ControlRGB',[1 0 0],'IsoRGB',[0 0 0],'MarkerCPT','o','MarkerCPTRGB',[1 0 0]);
			
				otherwise
					disp('Type not recognized')
            end
		end
                
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     function p = PolynomialOrder(obj)
         switch obj.type
             case 'curve'
                 p = obj.pu;
             case 'surf'
                 p = zeros(2,1);
                 p(1) = obj.pu;
                 p(2) = obj.pv;
             case 'volume'
                 p = zeros(3,1);
                 p(1) = obj.pu;
                 p(2) = obj.pv;
                 p(3) = obj.pw;
         end
     end
     function KnotCell = KnotVectorCell(obj)
         KnotCell = {};
         switch obj.type
             case 'curve'
                 KnotCell{1} = obj.U;
             case 'surf'
                 KnotCell{1} = obj.U;
                 KnotCell{2} = obj.V;
             case 'volume'
                 KnotCell{1} = obj.U;
                 KnotCell{2} = obj.V;
                 KnotCell{3} = obj.W;
         end   
     end

        function [PX, PY, PZ, weight] = UpdateCPTS(obj)
			switch obj.type
				case 'curve'
					PX = zeros(1,size(obj.Pw,2));
					PY = zeros(1,size(obj.Pw,2));
					PZ = zeros(1,size(obj.Pw,2));
					weight = zeros(1,size(obj.Pw,2));
				   for i=1:size(obj.Pw,2)
					   PX(i) = obj.Pw(i).x / obj.Pw(i).w;
					   PY(i) = obj.Pw(i).y / obj.Pw(i).w;
					   PZ(i) = obj.Pw(i).z / obj.Pw(i).w;
					   weight(i) = obj.Pw(i).w;
				   end
				   
				case 'surf'
					PX = zeros(size(obj.Pw,1),size(obj.Pw,2));
					PY = zeros(size(obj.Pw,1),size(obj.Pw,2));
					PZ = zeros(size(obj.Pw,1),size(obj.Pw,2));
					weight = zeros(size(obj.Pw,1),size(obj.Pw,2));

					for i=1:size(obj.Pw,1)
						for j=1:size(obj.Pw,2)
							PX(i,j) = obj.Pw(i,j).x / obj.Pw(i,j).w;
							PY(i,j) = obj.Pw(i,j).y / obj.Pw(i,j).w;
							PZ(i,j) = obj.Pw(i,j).z / obj.Pw(i,j).w;
							weight(i,j) = obj.Pw(i,j).w;
						end
					end

				case 'volume'
					PX = zeros(size(obj.Pw,1),size(obj.Pw,2),size(obj.Pw,3));
					PY = zeros(size(obj.Pw,1),size(obj.Pw,2),size(obj.Pw,3));
					PZ = zeros(size(obj.Pw,1),size(obj.Pw,2),size(obj.Pw,3));
					weight = zeros(size(obj.Pw,1),size(obj.Pw,2),size(obj.Pw,3));

					for i=1:size(obj.Pw,1)
						for j=1:size(obj.Pw,2)
							for k=1:size(obj.Pw,3)
								PX(i,j,k) = obj.Pw(i,j,k).x / obj.Pw(i,j,k).w;
								PY(i,j,k) = obj.Pw(i,j,k).y / obj.Pw(i,j,k).w;
								PZ(i,j,k) = obj.Pw(i,j,k).z / obj.Pw(i,j,k).w;
								weight(i,j,k) = obj.Pw(i,j,k).w;
							end
						end
					end
            end
        end
        
        
        
        function [PX, PY, PZ, weight] = UpdatePTS(obj,P)

			if iscell(P) ~= 1
				fprintf('P deve ser uma cell array\n');
			return
			end
			
			switch obj.type
				case 'curve'
					nu = size(P,2)-1;
					PX = zeros(1,nu+1);
					PY = zeros(1,nu+1);
					PZ = zeros(1,nu+1);
					weight = zeros(1,nu+1);


					for i=1:nu+1
						PX(i) = P{i}(1);
						PY(i) = P{i}(2);
						PZ(i) = P{i}(3);
						weight(i) = P{i}(4);
					end
					
				case 'surf'
					nu = size(P,1)-1;
					nv = size(P,2)-1;
					PX = zeros(nu+1,nv+1);
					PY = zeros(nu+1,nv+1);
					PZ = zeros(nu+1,nv+1);
					weight = zeros(nu+1,nv+1);

					for i=1:nu+1
						for j=1:nv+1
							PX(i,j) = P{i,j}(1);
							PY(i,j) = P{i,j}(2);
							PZ(i,j) = P{i,j}(3);
							weight(i,j) = P{i,j}(4);
						end
					end
					
				case 'volume'
					nu = size(P,1)-1;
					nv = size(P,2)-1;
					nw = size(P,3)-1;

					for i=1:nu+1
						for j=1:nv+1
							for k=1:nw+1
								PX(i,j,k) = P{i,j,k}(1);
								PY(i,j,k) = P{i,j,k}(2);
								PZ(i,j,k) = P{i,j,k}(3);
								weight(i,j,k) = P{i,j,k}(4);
							end
						end
					end
        	end
        end

        function [] = NormalizeKnotVector(obj)
            %very ugly, should create a separate function for a single knot vector later 
            switch obj.type
				case 'curve'
					if obj.U(1) < 0
						obj.U = obj.U - obj.U(1);
					elseif obj.U(1) > 0
						obj.U = obj.U + (0 - obj.U(1));
					end

					if obj.U(end) > 1
						obj.U = obj.U/obj.U(end);
					end
					
				case 'surf'
					if obj.U(1) < 0
						obj.U = obj.U - obj.U(1);
					elseif obj.U(1) > 0
						obj.U = obj.U + (0 - obj.U(1));
					end

					if obj.U(end) > 1
						obj.U = obj.U/obj.U(end);
					end
					%
					if obj.V(1) < 0
						obj.V = obj.V - obj.V(1);
					elseif obj.V(1) > 0
						obj.V = obj.V + (0 - obj.V(1));
					end

					if obj.V(end) > 1
						obj.V = obj.V/obj.V(end);
					end
                	%
					
				case 'volume'
                	%
					if obj.U(1) < 0
						obj.U = obj.U - obj.U(1);
					elseif obj.U(1) > 0
						obj.U = obj.U + (0 - obj.U(1));
					end

					if obj.U(end) > 1
						obj.U = obj.U/obj.U(end);
					end
					%
					if obj.V(1) < 0
						obj.V = obj.V - obj.V(1);
					elseif obj.V(1) > 0
						obj.V = obj.V + (0 - obj.V(1));
					end

					if obj.V(end) > 1
						obj.V = obj.V/obj.V(end);
					end
					%
					if obj.W(1) < 0
						obj.W = obj.W - obj.W(1);
					elseif obj.W(1) > 0
						obj.W = obj.W + (0 - obj.W(1));
					end

					if obj.W(end) > 1
						obj.W = obj.W/obj.W(end);
					end
            end 
        end
                
        function [] = set_PlotProp(obj,Color,LineSize,ControlColor,IsoColor,Marker,MarkerColor)
            
            switch nargin
                case 2
                    obj.PlotProp.RGB = Color;
                case 3
                    obj.PlotProp.RGB = Color;
                    obj.PlotProp.LineSize = LineSize;
                case 4
                    obj.PlotProp.RGB = Color;
                    obj.PlotProp.LineSize = LineSize;
                    obj.PlotProp.ControlRGB = ControlColor;
                case 5
                    obj.PlotProp.RGB = Color;
                    obj.PlotProp.LineSize = LineSize;
                    obj.PlotProp.ControlRGB = ControlColor;
                    obj.PlotProp.IsoRGB = IsoColor;
                case 6
                    obj.PlotProp.RGB = Color;
                    obj.PlotProp.LineSize = LineSize;
                    obj.PlotProp.ControlRGB = ControlColor;
                    obj.PlotProp.MarkerCPT = Marker;
                    obj.PlotProp.IsoRGB = IsoColor;
                    obj.PlotProp.MarkerCPT = Marker;
                case 7
                    obj.PlotProp.RGB = Color;
                    obj.PlotProp.LineSize = LineSize;
                    obj.PlotProp.ControlRGB = ControlColor;
                    obj.PlotProp.MarkerCPT = Marker;
                    obj.PlotProp.IsoRGB = IsoColor;
                    obj.PlotProp.MarkerCPT = Marker;
                    obj.PlotProp.MarkerCPTRGB = MarkerColor;
                otherwise
                    disp('Invalid input');
            end
            
        end
        
        function [] = reset_PlotProp(obj)
		
			switch obj.type
				case 'curve'
					obj.PlotProp = struct('RGB',[0 0 1],'LineSize',2,'ControlRGB',[0 0 0],'IsoRGB',[],'MarkerCPT','o','MarkerCPTRGB',[1 0 0]);
				case 'surf'
                	obj.PlotProp = struct('RGB',[.1 .9 .1],'LineSize',2,'ControlRGB',[1 0 0],'IsoRGB',[0 0 0],'MarkerCPT','o','MarkerCPTRGB',[1 0 0]);
				case 'volume'
                	obj.PlotProp = struct('RGB',[.1 .9 .1],'LineSize',2,'ControlRGB',[1 0 0],'IsoRGB',[0 0 0],'MarkerCPT','o','MarkerCPTRGB',[1 0 0]);
				otherwise
                	disp('Type not recognized');
			end
        end
                
                
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Geometry functions              
        function obj = KnotRefine(obj,X,dir)
            
			switch obj.type
				case 'curve'
					[obj.U, obj.Pw] = RefineKnotVectCurve(obj.nu,obj.pu,obj.U,obj.Pw,X,numel(X)-1);
					obj.nu = obj.nu + numel(X);
					[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);

				case 'surf'
					[obj.U, obj.V, obj.Pw] = RefineKnotVectSurface(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.Pw,X,numel(X)-1,dir);
					obj.nu = size(obj.Pw,1)-1;
					obj.nv = size(obj.Pw,2)-1;
					[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);
					
				case 'volume'
					if dir == 1
						[obj.U, obj.Pw] = RefineKnotSolid(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.nw,obj.pw,obj.W,obj.Pw,X,dir);
						obj.nu = size(obj.Pw,1)-1;
						[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);

					elseif dir == 2
						[obj.V, obj.Pw] = RefineKnotSolid(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.nw,obj.pw,obj.W,obj.Pw,X,dir);
						obj.nv = size(obj.Pw,2)-1;
						[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);

					elseif dir == 3
						[obj.W, obj.Pw] = RefineKnotSolid(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.nw,obj.pw,obj.W,obj.Pw,X,dir);
						obj.nw = size(obj.Pw,3)-1;
						[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);

					end
            end
        end
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
        function obj = KnotRemove(obj,knot,num,dir)
            
			switch obj.type
				case 'curve'
					r = FindSpanLinear(obj.nu,obj.pu,knot,obj.U);
					s = Mult(obj.nu,obj.pu,knot,obj.U);
					[~,obj.U,obj.Pw] = RemoveCurveKnot(obj.nu,obj.pu,obj.U,obj.Pw,knot,r,s,num);
					obj.nu = size(obj.Pw,2)-1;
					[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);
                
				case 'surf'
					if dir == 1
						r = FindSpanLinear(obj.nu,obj.pu,knot,obj.U);
						s = Mult(obj.nu,obj.pu,knot,obj.U);

					elseif dir == 2
						r = FindSpanLinear(obj.nv,obj.pv,knot,obj.V);
						s = Mult(obj.nv,obj.pv,knot,obj.V);
					end

					[~,obj.U,obj.V,obj.Pw] = RemoveSurfKnot(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.Pw,knot,dir,r,s,num);
					obj.nu = size(obj.Pw,1)-1;
					obj.nv = size(obj.Pw,2)-1;
					[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);
                
				case 'volume'
					if dir == 1
						[obj.U, obj.Pw] = RemoveKnotSolid(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.nw,obj.pw,obj.W,obj.Pw,knot,num,dir);
						obj.nu = size(obj.Pw,1)-1;
						[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);

					elseif dir == 2
						[obj.V, obj.Pw] = RemoveKnotSolid(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.nw,obj.pw,obj.W,obj.Pw,knot,num,dir);
						obj.nv = size(obj.Pw,2)-1;
						[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);

					elseif dir == 3
						[obj.W, obj.Pw] = RemoveKnotSolid(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.nw,obj.pw,obj.W,obj.Pw,knot,num,dir);
						obj.nw = size(obj.Pw,3)-1;
						[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);
					end
            end
        end
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
        function obj = DegreeElevate(obj,t,dir)
            
			switch obj.type
				case 'curve'
					[obj.nu,obj.U,obj.Pw] = DegreeElevateCurve(obj.nu,obj.pu,obj.U,obj.Pw,t);
					obj.pu = obj.pu + t;
					[obj.PX, obj.PY, obj.PZ, obj.weight] = obj.UpdateCPTS;
					
				case 'surf'
					[obj.nu,obj.U,obj.nv,obj.V,obj.Pw] = DegreeElevateSurface(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.Pw,dir,t);

					if dir == 1
						obj.pu = obj.pu + t;
					elseif dir == 2
						obj.pv = obj.pv + t;
					end

					[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);
					
				case 'volume'
					if dir == 1
						[obj.nu, obj.U, obj.Pw] = DegreeElevateSolid(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.nw,obj.pw,obj.W,obj.Pw,dir,t);
						obj.pu = obj.pu+t;
					elseif dir == 2
						[obj.nv, obj.V, obj.Pw] = DegreeElevateSolid(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.nw,obj.pw,obj.W,obj.Pw,dir,t);
						obj.pv = obj.pv+t;
					elseif dir == 3
						[obj.nw, obj.W, obj.Pw] = DegreeElevateSolid(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.nw,obj.pw,obj.W,obj.Pw,dir,t);
						obj.pw = obj.pw+t;
					end
					[obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);
            end
        end
                
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FEM connectivity data

        function [INN, IEN, nel, nen] = get_connectivity(obj)
			switch obj.type
				case 'curve'
					[INN, IEN, nel, nen] = GetConnectivity(obj.nu,obj.pu);
				case 'surf'
					[INN, IEN, nel, nen] = GetConnectivity(obj.nu,obj.pu,obj.nv,obj.pv);
				case 'volume'
					[INN, IEN, nel, nen] = GetConnectivity(obj.nu,obj.pu,obj.nv,obj.pv,obj.nw,obj.pw);
				otherwise
					disp('Invalid type')
            end
        end
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plots
        function [] = plot_geo(obj,render,cpoints,isolines,Urange,Vrange,Wrange)
            
            if strcmp(obj.type,'curve')
                
                switch nargin
                    case 1
                        render = 'medium';
                        cpoints = 0;
                        isolines = 0;
                        Urange = [obj.U(1) obj.U(end)];
                    case 2
                        cpoints = 0;
                        isolines = 0;
                        Urange = [obj.U(1) obj.U(end)];
                    case 3
                        isolines = 0;
                        Urange = [obj.U(1) obj.U(end)];
                    case 4
                        Urange = [obj.U(1) obj.U(end)];
                end
                
                PlotRatCurve(obj.nu,obj.pu,obj.U,obj.Pw,Urange,cpoints,isolines,obj.PlotProp);
                    
            
            elseif strcmp(obj.type,'surf')
                
                switch nargin
                    case 1
                        render = 'medium';
                        cpoints = 0;
                        isolines = 0;
                        Urange = [obj.U(1) obj.U(end)];
                        Vrange = [obj.V(1) obj.V(end)];
                    case 2
                        cpoints = 0;
                        isolines = 0;
                        Urange = [obj.U(1) obj.U(end)];
                        Vrange = [obj.V(1) obj.V(end)];
                    case 3
                        isolines = 0;
                        Urange = [obj.U(1) obj.U(end)];
                        Vrange = [obj.V(1) obj.V(end)];
                    case 4
                        Urange = [obj.U(1) obj.U(end)];
                        Vrange = [obj.V(1) obj.V(end)];
                    case 5
                        Vrange = [obj.V(1) obj.V(end)];
                end
                        
                PlotSurf(obj,render,Urange,Vrange);
                hold on;
                
                if isolines
                    Unique = unique(obj.U);
                    Vnique = unique(obj.V);
                    for i=1:numel(Unique)
                        PlotSurfPatches_1(obj.PX,obj.PY,obj.PZ,obj.weight,obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,render,Unique(i),obj.PlotProp);
                    end
    
                    for j=1:numel(Vnique)
                        PlotSurfPatches_2(obj.PX,obj.PY,obj.PZ,obj.weight,obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,render,Vnique(j),obj.PlotProp);
                    end
                end
                
                if cpoints
                    plot3(obj.PX,obj.PY,obj.PZ,'color',obj.PlotProp.ControlRGB,'LineWidth',obj.PlotProp.LineSize);
                    plot3((obj.PX)',(obj.PY)',(obj.PZ)','color',obj.PlotProp.ControlRGB,'LineWidth',obj.PlotProp.LineSize);
                    plot3(obj.PX,obj.PY,obj.PZ,obj.PlotProp.MarkerCPT,'MarkerEdgeColor',obj.PlotProp.MarkerCPTRGB, 'MarkerFaceColor',obj.PlotProp.MarkerCPTRGB);
                end
                light;

                
            elseif strcmp(obj.type,'volume')
                
                switch nargin
                    case 1  
                        render = 'medium';
                        cpoints = 0;
                        isolines = 0;
                        Urange = [obj.U(1) obj.U(end)];
                        Vrange = [obj.V(1) obj.V(end)];
                        Wrange = [obj.W(1) obj.W(end)];
                    case 2
                        cpoints = 0;
                        isolines = 0;
                        Urange = [obj.U(1) obj.U(end)];
                        Vrange = [obj.V(1) obj.V(end)];
                        Wrange = [obj.W(1) obj.W(end)];
                    case 3
                        isolines = 0;
                        Urange = [obj.U(1) obj.U(end)];
                        Vrange = [obj.V(1) obj.V(end)];
                        Wrange = [obj.W(1) obj.W(end)];
                    case 4
                        Urange = [obj.U(1) obj.U(end)];
                        Vrange = [obj.V(1) obj.V(end)];
                        Wrange = [obj.W(1) obj.W(end)];
                    case 5
                        Vrange = [obj.V(1) obj.V(end)];
                        Wrange = [obj.W(1) obj.W(end)];
                    case 6
                        Wrange = [obj.W(1) obj.W(end)];
                end
                
                if nargin > 4
                    PlotRatVolRange(obj,Urange,Vrange,Wrange,render,cpoints,isolines)
                else
                    PlotRatVol(obj,render,cpoints,isolines);
                end
                
                xlabel('x [m]','FontWeight','bold','FontSize',23)
                ylabel('y [m]','FontWeight','bold','FontSize',23)
                zlabel('z [m]','FontWeight','bold','FontSize',23)
                set(gca,'FontSize',23)
                set(gca,'FontWeight','bold')
            end
            
        end
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
        function [] = plot_basis(obj,dir)
            
			switch obj.type
				case 'curve'
                	PlotAllBasisFuns(obj.nu,obj.pu,obj.U);
				case 'surf'
					switch dir
						case 1
							PlotAllBasisFuns(obj.nu,obj.pu,obj.U);
						case 2
							PlotAllBasisFuns(obj.nv,obj.pv,obj.V);
						otherwise
							disp('Invalid direction')
					end
                
				case 'volume'
					switch dir
						case 1
							PlotAllBasisFuns(obj.nu,obj.pu,obj.U);
						case 2
							PlotAllBasisFuns(obj.nv,obj.pv,obj.V);
						case 3
							PlotAllBasisFuns(obj.nw,obj.pw,obj.W);
						otherwise
							disp('Invalid direction')
					end 
            end
        end
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        function point = eval_point(obj,knot,knot2,knot3)
            
			switch obj.type
				case 'curve'
                	point = CCurvePoint2(obj.nu,obj.pu,obj.U,obj.Pw,knot);
            	case 'surf'
                	S = SurfacePointRAT3(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.Pw,knot,knot2);
                	point = CPOINT(S.x,S.y,S.z,1,1);
            	case 'volume'
                	point = VolumePoint(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.nw,obj.pw,obj.W,obj.Pw,knot,knot2,knot3);
			end
        end
		
		function dR = eval_derivative(obj,u,v,P)
		
			switch obj.type
				case 'curve'
					disp('yet to be implemented')
				case 'surf'
					dR = CalculateDers2D(obj,u,v,P);
				case 'volume'
					disp('yet to be implemented')
				otherwise
					disp('Invalid type')
			end
		end
					

        function Nip = eval_basis(obj,dir,i,knot)
				switch dir
					case 1
						Nip = OneBasisFun(obj.pu,obj.nu +obj.pu +1,obj.U,i,knot);
					case 2
						Nip = OneBasisFun(obj.pv,obj.nv+obj.pv+1,obj.V,i,knot);
					case 3
						Nip = OneBasisFun(obj.pw,obj.nw+obj.pw+1,obj.W,i,knot);
				end
		end
		
        function P = get_point_cell(obj)
			
			switch obj.type
				case 'volume'
					P = cell(size(obj.Pw,1),size(obj.Pw,2),size(obj.Pw,3));
               
					nu = size(obj.Pw,1); nv = size(obj.Pw,2); nw = size(obj.Pw,3);
               
					for i=1:nu
						for j=1:nv
							for k=1:nw
								P{i,j,k}(1) = obj.PX(i,j,k);
								P{i,j,k}(2) = obj.PY(i,j,k);
								P{i,j,k}(3) = obj.PZ(i,j,k);
								P{i,j,k}(4) = obj.weight(i,j,k);
							end
						end
					end
					
				case 'curve'
					P = cell(length(obj.Pw),1);
					nu = length(obj.Pw);

					for i=1:nu
						P{i}(1) = obj.PX(i);
						P{i}(2) = obj.PY(i);
						P{i}(3) = obj.PZ(i);
						P{i}(4) = obj.weight(i);
					end
					
				case 'surf'

					P = cell(size(obj.Pw,1),size(obj.Pw,2));
					nu = size(obj.Pw,1); nv = size(obj.Pw,2);
               
					for i=1:nu
						for j=1:nv
							P{i,j}(1) = obj.PX(i,j);
							P{i,j}(2) = obj.PY(i,j);
							P{i,j}(3) = obj.PZ(i,j);
							P{i,j}(4) = obj.weight(i,j);
						end
					end
			end 
        end		
 
    end
end

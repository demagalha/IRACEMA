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
            
            if strcmp(varargin{1},'curve')
                
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
            end
            
            if strcmp(varargin{1},'surf')

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
            end
            
            if strcmp(varargin{1},'volume')
                
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
            end
        end
                
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
        function [PX, PY, PZ, weight] = UpdateCPTS(obj)
            
            if strcmp(obj.type,'curve') == 1
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
            
            
            elseif strcmp(obj.type,'surf') == 1
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
                
            elseif strcmp(obj.type,'volume') == 1
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

    
        if strcmp(obj.type,'curve') == 1
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
        
        elseif strcmp(obj.type,'surf') == 1
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
        
        elseif strcmp(obj.type,'volume') == 1
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
            
            if strcmp(obj.type,'curve')
                
                if obj.U(1) < 0
                    obj.U = obj.U - obj.U(1);
                elseif obj.U(1) > 0
                    obj.U = obj.U + (0 - obj.U(1));
                end
                
                if obj.U(end) > 1
                    obj.U = obj.U/obj.U(end);
                end
                
            elseif strcmp(obj.type,'surf')
                
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
            elseif strcmp(obj.type,'volume')
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
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
        function obj = KnotRefine(obj,X,dir)
            
            if strcmp(obj.type,'curve') == 1
                [obj.U, obj.Pw] = RefineKnotVectCurve(obj.nu,obj.pu,obj.U,obj.Pw,X,numel(X)-1);
                obj.nu = obj.nu + numel(X);
                [obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);
                
            
            elseif strcmp(obj.type,'surf') == 1
                [obj.U, obj.V, obj.Pw] = RefineKnotVectSurface(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.Pw,X,numel(X)-1,dir);
                obj.nu = size(obj.Pw,1)-1;
                obj.nv = size(obj.Pw,2)-1;
                [obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);
                
                
            elseif strcmp(obj.type,'volume') == 1
                
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
            
            if strcmp(obj.type,'curve') == 1
                r = FindSpanLinear(obj.nu,obj.pu,knot,obj.U);
                s = Mult(obj.nu,obj.pu,knot,obj.U);
                [~,obj.U,obj.Pw] = RemoveCurveKnot(obj.nu,obj.pu,obj.U,obj.Pw,knot,r,s,num);
                obj.nu = size(obj.Pw,2)-1;
                [obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);
                
            elseif strcmp(obj.type, 'surf') == 1
                
                if dir == 1
                 r = FindSpanLinear(obj.nu,obj.pu,knot,obj.U);
                 s = Mult(obj.nu,obj.pu,knot,obj.U);
                 
                elseif dir ==2
                    r = FindSpanLinear(obj.nv,obj.pv,knot,obj.V);
                    s = Mult(obj.nv,obj.pv,knot,obj.V);
                end
            
                [~,obj.U,obj.V,obj.Pw] = RemoveSurfKnot(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.Pw,knot,dir,r,s,num);
                obj.nu = size(obj.Pw,1)-1;
                obj.nv = size(obj.Pw,2)-1;
                [obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);
                
            elseif strcmp(obj.type,'volume') == 1
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
            
            if strcmp(obj.type,'curve') == 1
                [obj.nu,obj.U,obj.Pw] = DegreeElevateCurve(obj.nu,obj.pu,obj.U,obj.Pw,t);
                obj.pu = obj.pu + t;
                [obj.PX, obj.PY, obj.PZ, obj.weight] = obj.UpdateCPTS;
                
            elseif strcmp(obj.type,'surf') == 1
                [obj.nu,obj.U,obj.nv,obj.V,obj.Pw] = DegreeElevateSurface(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.Pw,dir,t);
				
				if dir == 1
				obj.pu = obj.pu + t;
				elseif dir == 2
				obj.pv = obj.pv + t;
				end
				
                [obj.PX, obj.PY, obj.PZ, obj.weight] = UpdateCPTS(obj);
                
            elseif strcmp(obj.type,'volume') == 1
                
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
                
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        function [INN, IEN, nel, nen] = get_connectivity(obj)
            if strcmp(obj.type,'curve') == 1
                [INN, IEN, nel, nen] = GetConnectivity(obj.nu,obj.pu);
            elseif strcmp(obj.type,'surf') == 1
                [INN, IEN, nel, nen] = GetConnectivity(obj.nu,obj.pu,obj.nv,obj.pv);
            else
                [INN, IEN, nel, nen] = GetConnectivity(obj.nu,obj.pu,obj.nv,obj.pv,obj.nw,obj.pw);
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = plot_geo(obj,render,cpoints,isolines)
            
            if strcmp(obj.type,'curve')
                if nargin <= 3
                    cpoints = 0;
                    isolines = 0;
                end
                
                PlotRatCurve(obj.nu,obj.pu,obj.U,obj.Pw,cpoints,isolines,obj.PlotProp);
                    
            
            elseif strcmp(obj.type,'surf')
                
                if (nargin <= 3)
                    render = 'medium';
                    cpoints = 0;
                    isolines = 0;
                    
                end
                PlotSurf(obj.PX,obj.PY,obj.PZ,obj.weight,obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,render,obj.PlotProp.RGB);
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
                
                if (nargin == 1)
                    render = 'medium';
                    cpoints = 0;
                    isolines = 0;
                    
                end
                PlotRatVol(obj.PX,obj.PY,obj.PZ,obj.weight,obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.nw,obj.pw,obj.W,render,cpoints,isolines,obj.PlotProp);

                xlabel('x [m]','FontWeight','bold','FontSize',23)
                ylabel('y [m]','FontWeight','bold','FontSize',23)
                zlabel('z [m]','FontWeight','bold','FontSize',23)
                set(gca,'FontSize',23)
                set(gca,'FontWeight','bold')
            end
            
        end
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
        function [] = plot_basis(obj,dir)
            
            if strcmp(obj.type,'curve') == 1
                PlotAllBasisFuns(obj.nu,obj.pu,obj.U);
                
            elseif strcmp(obj.type,'surf') == 1
                if dir == 1
                    PlotAllBasisFuns(obj.nu,obj.pu,obj.U);
                else
                    PlotAllBasisFuns(obj.nv,obj.pv,obj.V);
                end
                
            elseif strcmp(obj.type,'volume') == 1
                if dir == 1
                    PlotAllBasisFuns(obj.nu,obj.pu,obj.U);
                elseif dir == 2
                    PlotAllBasisFuns(obj.nv,obj.pv,obj.V);
                else
                    PlotAllBasisFuns(obj.nw,obj.pw,obj.W);
                end
                
            end
        end
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        function point = eval_point(obj,knot,knot2,knot3)
            
            if strcmp(obj.type,'curve')
                point = CCurvePoint2(obj.nu,obj.pu,obj.U,obj.Pw,knot);
                
            elseif strcmp(obj.type,'surf')
                S = SurfacePointRAT3(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.Pw,knot,knot2);
                point = CPOINT(S.x,S.y,S.z,1,1);
                
            elseif strcmp(obj.type,'volume')
                point = VolumePoint(obj.nu,obj.pu,obj.U,obj.nv,obj.pv,obj.V,obj.nw,obj.pw,obj.W,obj.Pw,knot,knot2,knot3);
            end
            
        end

        function Nip = eval_basis(obj,dir,i,knot)
               if dir == 1
                Nip = OneBasisFun(obj.pu,obj.nu +obj.pu +1,obj.U,i,knot);
               elseif dir == 2
                   Nip = OneBasisFun(obj.pv,obj.nv+obj.pv+1,obj.V,i,knot);
               elseif dir == 3
                   Nip = OneBasisFun(obj.pw,obj.nw+obj.pw+1,obj.W,i,knot);
               end
        end
		
        function P = get_point_cell(obj)
           
            if strcmp('volume',obj.type) == 1
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
           
            elseif strcmp('curve',obj.type) == 1
                P = cell(length(obj.Pw),1);
                nu = length(obj.Pw);
               
                for i=1:nu
                    P{i}(1) = obj.PX(i);
                    P{i}(2) = obj.PY(i);
                    P{i}(3) = obj.PZ(i);
                    P{i}(4) = obj.weight(i);
                end
           
            elseif strcmp('surf',obj.type) == 1
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
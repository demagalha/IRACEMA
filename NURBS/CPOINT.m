classdef CPOINT
    properties
        x;
        y;
        z;
        w;
    end
    methods
        function obj = CPOINT(x,y,z,w,point)
            if point == 1
                obj.x = x;
                obj.y = y;
                obj.z = z;
                obj.w = w;
            else
                
            obj.x = x .* w;
            obj.y = y .* w;
            obj.z = z .*w;
            obj.w = w;
            end
        end
        
        function D = Distance4D(P1,P2)
            D = (P2.x - P1.x)^2 + (P2.y - P1.y)^2 + (P2.z - P1.z)^2 + (P2.w - P1.w)^2;
            D = sqrt(D);
        end
        
        function r = plus(a,b)
            if isobject(a) && isobject(b)
                r = CPOINT(a.x+b.x,a.y+b.y,a.z+b.z,a.w+b.w,1);
            end
            if isobject(a) && (isobject(b) ~= 1)
                r = CPOINT(a.x+b,a.y+b,a.z+b,a.w+b,1);
            end
            %if isobject(a) ~=1 && isobject(b)
             %   r = CPOINT(a+b.x,a+b.y,a+b.z,a+b.w,1);
            %end
        end
        function s = minus(a,b)
            if isobject(a) && isobject(b)
                s = CPOINT(a.x-b.x,a.y-b.y,a.z-b.z,a.w-b.w,1);
            end
            if isobject(a) && (isobject(b) ~= 1)
                s = CPOINT(a.x-b,a.y-b,a.z-b,a.w-b,1);
            end
            %if isobject(a) ~=1 && isobject(b)
             %   s = CPOINT(a-b.x,a-b.y,a-b.z,a-b.w,1);
            %end
        end
        function t = mtimes(a,b)
            %%%%elemento por elemento
            if isobject(a) && isobject(b)
                t = CPOINT(a.x*b.x,a.y*b.y,a.z*b.z,a.w*b.w,1);
            end
            %%%%%
            if isobject(a) && (isobject(b) ~= 1)
                t = CPOINT(a.x*b,a.y*b,a.z*b,a.w*b,1);
            end
            if isobject(a) ~=1 && isobject(b)
                t = CPOINT(a*b.x,a*b.y,a*b.z,a*b.w,1);
            end
        end
        function u  = mrdivide(a,b)
            
            if isobject(a) && isobject(b)
                u = CPOINT(a.x/b.x,a.y/b.y,a.z/b.z,a.w/b.w,1);
            end
            %%%%%
            if isobject(a) && (isobject(b) ~= 1)
                u = CPOINT(a.x/b,a.y/b,a.z/b,a.w/b,1);
            end
            %if isobject(a) ~=1 && isobject(b)
             %   u = CPOINT(a/b.x,a/b.y,a/b.z,a/b.w,1);
            %end
        end

    end
end
            
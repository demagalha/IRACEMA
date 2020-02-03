function [CGPoints, CGPointsMinus] = UnivariateCauchyGalerkin(KnotVector, p)
check = rem(p,2);
tmp = unique(KnotVector);
    if logical(check)
        if p==3
            x = 1/sqrt(3);
        elseif p==5
            x = sqrt(225 -30*sqrt(30))/15;
        elseif p==7
            x = 0.5049185675126533;
        end
            for i=1:numel(tmp)-1
                b = tmp(i+1);
                a = tmp(i);
                CGPoints(i) = ((b-a)/2)*x + (b+a)/2;
                CGPointsMinus(i) = ((b-a)/2)*(-x) + (b+a)/2;
            end
    else
        x = 0.5;
        for i=1:numel(tmp)-1
            b = tmp(i+1);
            a = tmp(i);
            CGPoints(i) = (b+a)/2;
        end
        CGPointsMinus = tmp;
    end
        
end
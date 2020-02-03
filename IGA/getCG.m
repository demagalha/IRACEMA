function [u, wu] = getCG(p)
check = rem(p,2);
if logical(check)
    if (p == 3) || (p == 1)
        u(1) = 1/sqrt(3);
        u(2) = -u(1);
    elseif p == 5
        u(1) = sqrt(225 -30*sqrt(30))/15;
        u(2) = -u(1);
    elseif p==7
        u(1) = 0.5049185675126533;
        u(2) = -u(1);
    end
    wu = ones(2,1);
else
    u(1) = -1;
    u(2) = 0;
    u(3) = 1;
    wu = ones(3,1);
end
end


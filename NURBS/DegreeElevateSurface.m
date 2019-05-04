function [nh,Uh,mh,Vh,Qw] = DegreeElevateSurface(n,p,U,m,q,V,Pw,dir,t)

if dir == DIR.U_DIRECTION
    s = FindS(U,p);
    Qw(1:n+t*(1+s)+1,1:m+1) = CPOINT(0,0,0,0,1);
    
    for row=0:m
        [nh,Uh,Qw(:,row+1)] = DegreeElevateCurve(n,p,U,Pw(:,row+1),t);
    end
    Vh = V;
    mh = m;
end

if dir == DIR.V_DIRECTION
    s = Find(V,q);
    Qw(1:n+1,1:m+t*(1+s)+1) = CPOINT(0,0,0,0,1);
    
    for row=0:n
        [mh,Vh,Qw(row+1,:)] = DegreeElevateCurve(n,p,U,Pw(row+1,:),t);
    end
    Uh = U;
    nh = n;
end

end
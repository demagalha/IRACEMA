function [Ubar, Vbar, Qw] = RefineKnotVectSurface(n,p,U,m,q,V,Pw,X,r,dir)

if dir == 1
    
    Qw = Pw;
    for i=0:r
    k = FindSpanLinear(n,p,X(i+1),U);
    s = Mult(n,p,X(i+1),U);    
    [n,U,m,V,Qw] = SurfaceKnotIns(n,p,U,m,q,V,Qw,dir,X(i+1),k,s,1);
    end
end

if dir == 2
    
    Qw = Pw;
    for i=0:r
        k = FindSpanLinear(m,q,X(i+1),V);
        s = Mult(m,q,X(i+1),V);
        [n,U,m,V,Qw] = SurfaceKnotIns(n,p,U,m,q,V,Qw,dir,X(i+1),k,s,1);
    end
end

Ubar = U;
Vbar = V;


end
function [Unew, Qw] = RefineKnotSolid(nu,pu,U,nv,pv,V,nw,pw,W,Pw,X,dir)

%%%knot insetion to a NURBS solid %%% tirei r antes do dir (inutil)

if dir == 1 
Qw(1:size(Pw,1)+length(X),1:size(Pw,2),1:size(Pw,3)) = CPOINT(0,0,0,0,0);

P(1:size(Pw,1)) = CPOINT(0,0,0,0,1);

    
    for j=1:size(Pw,2)
        for k=1:size(Pw,3)
            for i=1:size(Pw,1)
                P(i) = CPOINT(Pw(i,j,k).x,Pw(i,j,k).y,Pw(i,j,k).z,Pw(i,j,k).w,1);
            end
            r = length(X)-1;
            [Ubar, Q] = RefineKnotVectCurve(size(P,2)-1,pu,U,P,X,r);
         
            for i=1:size(Pw,1)+length(X)
             Qw(i,j,k) = Q(i);
            end
        end
    end

Unew = Ubar;
end

if dir == 2 %REFINE IN THE V DIRECTION
    Qw(1:size(Pw,1),1:size(Pw,2)+length(X),1:size(Pw,3)) = CPOINT(0,0,0,0,0);

    P(1:size(Pw,2)) = CPOINT(0,0,0,0,1);

    for i=1:size(Pw,1)
        for k=1:size(Pw,3)
            for j=1:size(Pw,2)
                P(j) = CPOINT(Pw(i,j,k).x,Pw(i,j,k).y,Pw(i,j,k).z,Pw(i,j,k).w,1);
            end
            r = length(X)-1;
            [Vbar, Q] = RefineKnotVectCurve(size(P,2)-1,pv,V,P,X,r);
         
            for j=1:size(Pw,2)+length(X)
             Qw(i,j,k) = Q(j);
            end
        end
    end
Unew = Vbar;
end

if dir == 3 %%%W DIRECTION
    Qw(1:size(Pw,1),1:size(Pw,2),1:size(Pw,3)+length(X)) = CPOINT(0,0,0,0,0);
    
    P(1:size(Pw,3)) = CPOINT(0,0,0,0,1);
    
    for i=1:size(Pw,1)
        for j=1:size(Pw,2)
            for k=1:size(Pw,3)
                P(k) = CPOINT(Pw(i,j,k).x,Pw(i,j,k).y,Pw(i,j,k).z,Pw(i,j,k).w,1);
            end
            r = length(X)-1;
            [Wbar, Q] = RefineKnotVectCurve(size(P,2)-1,pw,W,P,X,r);
         
            for k=1:size(Pw,3)+length(X)
             Qw(i,j,k) = Q(k);
            end
        end
    end
Unew = Wbar;
end

end

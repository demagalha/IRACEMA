function [Unew, Qw] = RemoveKnotSolid(nu,pu,U,nv,pv,V,nw,pw,W,Pw,knot,num,dir)

if dir == 1
    
Qw(1:size(Pw,1)-num,1:size(Pw,2),1:size(Pw,3)) = CPOINT(0,0,0,0,0);    
P(1:size(Pw,1)) = CPOINT(0,0,0,0,1);
    %%%adicionar um check para ver se o knot é removível (para direção U no
    %%%caso size(Pw,2)*size(Pw,3) vezes, se soma t ~= size(Pw,2)*size(P2,3)
    %%%... sai da função, já que não é possível remover, ou exit function
    %%%assim que não der pra remover um knot (dá no mesmo)
    
     for j=1:size(Pw,2)
        for k=1:size(Pw,3)
            for i=1:size(Pw,1)
                P(i) = CPOINT(Pw(i,j,k).x,Pw(i,j,k).y,Pw(i,j,k).z,Pw(i,j,k).w,1);
            end
            r = FindSpanLinear(size(P,2)-1,pu,knot,U);
            s = Mult(size(P,2)-1,pu,knot,U);
            [t,Ubar,Q] = RemoveCurveKnot(size(P,2)-1,pu,U,P,knot,r,s,num);
         
            for i=1:size(Pw,1)-num
             Qw(i,j,k) = Q(i);
            end
        end
     end
    Unew = Ubar;
end

if dir == 2

Qw(1:size(Pw,1),1:size(Pw,2)-num,1:size(Pw,3)) = CPOINT(0,0,0,0,0);    
P(1:size(Pw,2)) = CPOINT(0,0,0,0,1);

 for i=1:size(Pw,1)
        for k=1:size(Pw,3)
            for j=1:size(Pw,2)
                P(j) = CPOINT(Pw(i,j,k).x,Pw(i,j,k).y,Pw(i,j,k).z,Pw(i,j,k).w,1);
            end
            r = FindSpanLinear(size(P,2)-1,pv,knot,V);
            s = Mult(size(P,2)-1,pv,knot,V);
            [t,Vbar,Q] = RemoveCurveKnot(size(P,2)-1,pv,V,P,knot,r,s,num);
         
            for j=1:size(Pw,2)-num
             Qw(i,j,k) = Q(j);
            end
        end
     end
    Unew = Vbar;
end

if dir == 3
    
    Qw(1:size(Pw,1),1:size(Pw,2),1:size(Pw,3)-num) = CPOINT(0,0,0,0,0);    
    P(1:size(Pw,3)) = CPOINT(0,0,0,0,1);
    
     for i=1:size(Pw,1)
        for j=1:size(Pw,2)
            for k=1:size(Pw,3)
                P(k) = CPOINT(Pw(i,j,k).x,Pw(i,j,k).y,Pw(i,j,k).z,Pw(i,j,k).w,1);
            end
            r = FindSpanLinear(size(P,2)-1,pw,knot,W);
            s = Mult(size(P,2)-1,pw,knot,W);
            [t,Wbar,Q] = RemoveCurveKnot(size(P,2)-1,pw,W,P,knot,r,s,num);
         
            for k=1:size(Pw,3)-num
             Qw(i,j,k) = Q(k);
            end
        end
     end
    Unew = Wbar;
end

end
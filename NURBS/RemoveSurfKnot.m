function [t,Unew,Vnew,Pww] = RemoveSurfKnot(n,p,U,m,q,V,Pw,u,dir,r,s,num)


%%%%%%%%%%%%%%%%%%%%%como se um nó não for removível já não da pra remover,
%%%%%%%%%%%%%%%%%%%%%colocar um while(flag ==1 && t(row+1) == 1), <- evita
%%%%%%%%%%%%%%%%%%%%%ir até o fim

if dir == DIR.U_DIRECTION || dir == 1
    %Pw(:,1)
    t = zeros(1,m+1);
    
    for row=0:m
        %[Pw(:,row+1).x]
        [t(row+1),~,~] = RemoveCurveKnot(n,p,U,Pw(:,row+1),u,r,s,num);
    end
  
    if sum(t) ~= m+1
        fprintf('Cannot remove knot\n exiting function')
    end
        
end
           
if dir == DIR.V_DIRECTION || dir == 2
    %Pw(1,:)
    t = zeros(1,n+1);
    
    %não é row (é columm na direção V) per se, row é pro U_DIRECTION mas anyway
    for row=0:n
        %[Pw(row+1,:).x]
        [t(row+1),~,~] = RemoveCurveKnot(m,q,V,Pw(row+1,:),u,r,s,num);
    end
    
    if sum(t) ~= n+1
        fprintf('Cannot remove knot\n exiting function')
        return;
    end
end

%dir


if dir == DIR.U_DIRECTION || dir == 1
    Pww(1:n+1-num,1:m+1) = CPOINT(0,0,0,1,1);
    
    for row=0:m
        [~,Unew,Pww(:,row+1)] = RemoveCurveKnot(m,q,V,Pw(:,row+1),u,r,s,num);
    end
    Vnew = V;
end

if dir == DIR.V_DIRECTION || dir == 2
    
    Pww(1:n+1,1:m+1-num) = CPOINT(0,0,0,1,1);
    for row=0:n
        %[Pw(row+1,:).x]
        
        %fprintf('teste')
        [~,Vnew,Pww(row+1,:)] = RemoveCurveKnot(m,q,V,Pw(row+1,:),u,r,s,num);
    end
    Unew = U;
end



end
        
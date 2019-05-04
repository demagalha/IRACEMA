function [t,U,Pw] = RemoveCurveKnot(n,p,U,Pw,u,r,s,num)

TOL = .75;

m = n+p+1;
ord = p+1;
fout = floor((2*r-s-p)/2);
last = r-s;
first = r-p;

temp(1:2*p+1) = CPOINT(0,0,0,0,1);
t = 0;
while t<num
    off = first-1;
    temp(1) = Pw(off+1);
    temp(last+1-off+1) = Pw(last+2);
    i = first;
    j = last;
    ii = 1;
    jj = last-off;
    remflag = 0;
    while j-i > t
        alfi = (u-U(i+1))/(U(i+ord+t+1)-U(i+1));
        alfj = (u-U(j-t+1))/(U(j+ord+1)-U(j-t+1));
        temp(ii+1) = (Pw(i+1)-(1-alfi)*temp(ii))/alfi;
        temp(jj+1) = (Pw(j+1)-alfj*temp(jj+2))/(1-alfj);
        i = i +1;
        ii = ii+1;
        j = j-1;
        jj = jj-1;
    end
    
    if j-i < t
        if (Distance4D(temp(ii),temp(jj+2)) <= TOL)
            %fprintf('ola')
            remflag = 1;
        end
    else
        alfi = (u-U(i+1))/(U(i+ord+t+1)-U(i+1));
        if (Distance4D(Pw(i+1),alfi*temp(ii+t+2)+(1-alfi)*temp(ii)) <= TOL)
            remflag = 1;
            %fprintf('ola2')
        end
    end
    
    if remflag == 0 %não consegue remover mais
        break       %saindo do loop
    else
        %removido com sucesso, salvar cpts
        i = first;
        j = last;
        while j-i > t
            Pw(i+1) = temp(i-off+1);
            Pw(j+1) = temp(j-off+1);
            i = i + 1;
            j = j - 1;
        end
    end
    
    first = first - 1;
    last = last +1;
    t = t+1;
end %end of for

if t == 0
    return;
end

for k =r+1:m
    U(k-t+1) = U(k+1);
end
    j = fout;
    i = j;
    
    for k=1:t-1
        if mod(k,2) == 1
            i = i +1;
        else
            j = j-1;
        end
    end
    
    for k=i+1:n
        Pw(j+1) = Pw(k+1);
        j = j+1;
    end
U = U(1:end-t);
Pw = Pw(1:end-t);
    return;
end
        
        
        
        
        
        
        
        
        
function s = FindS(U,p)


if length(U) == 2*(p+1)
    s = 0;
    return;
end

Uaux = U(p+2:end-(p+1));

tam = length(Uaux);

if tam == 1
    s =1;
    return;
end

s = 1;
for i=2:tam
    if(Uaux(i) ~= Uaux(i-1))
        s = s + 1;
    end
end
end
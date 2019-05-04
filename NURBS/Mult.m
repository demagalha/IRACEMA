function s = Mult(n,p,u,U)


s = 0;

if(u == U(n+2))
    s = p+1;
    return
end

if(u == U(p+1))
    s = p+1;
    return
end

for i=p+2:n+1
    if(u == U(i))
        s = s+1;
    end
end

end
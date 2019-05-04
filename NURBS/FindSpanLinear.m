function KnotSpanLinear = FindSpanLinear(n, p, u, U)

%implementação para determinar o índice do knot span por uma busca linear
%entrada: n, p, u, U 
%U(knot vector), u valor arbitrário, p ordem e n = m-p-1
%saída índice do knot span

%m = length(U);


    if u == U(n+2)
        KnotSpanLinear = n;
        return
    end
    
    for j=1:n+1
        if(u >= U(j) && u < U(j+1))
            KnotSpanLinear = j -1;
            return
        end
    end
end

    
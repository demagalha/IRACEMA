function [knots_1, knots_2] = KnotToInsert(U_1,U_2)

a = U_1(1);
b = U_2(1);
p_1 = 0;
p_2 = 0;

knots_1 = 0;
knots_2 = 0;
for i=1:numel(U_1)
    if U_1(i) ~= a
        break;
    end
    p_1 = p_1 + 1;
end

for i=1:numel(U_2)
    if U_2(i) ~= b
        break;
    end
    p_2 = p_2 + 1;
end

if p_1 ~= p_2
    fprintf('Knot vectors must be of same degree');
    return;
end
V_1 = U_1;
V_2 = U_2;
p_1 = p_1-1;
p_2 = p_2-1;

%if numel(U_1) > numel(U_2)
    for i = p_1+1+1:numel(U_1)-(p_1+1)
        for j = p_2+1+1:numel(U_2)-(p_2+1)
            
            if V_1(i) == V_2(j)
                V_1(i) = -1;
                V_2(j) = -1;
                break;

            end
        end
    end
    
    for i = p_2+1+1:numel(U_2)-(p_2+1)
        for j = p_1+1+1:numel(U_1)-(p_1+1)
            
            if U_2(i) == U_1(j)
                U_2(i) = -1;
                U_1(j) = -1;
            end
        end
    end
%end



a = 1;
b = 1;
for i = p_1+1+1:numel(V_1)-(p_1+1)
    if V_1(i) ~= -1
        knots_1(a) = V_1(i);
        a = a+1;
    end
    
end

for i = p_2+1+1:numel(U_2)-(p_2+1)
    if U_2(i) ~= -1
        knots_2(b) = U_2(i);
        b = b+1;
    end

end

                
                
       
end    
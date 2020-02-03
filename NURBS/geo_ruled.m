function ruled = geo_ruled(curve_1,curve_2,invert_1,invert_2)


p_1 = curve_1.pu;
p_2 = curve_2.pu;

p = max(p_1,p_2);




if p_1 > p_2
    curve_2.DegreeElevate(p_1 - p_2,1)
elseif p_2 > p_1
    curve_1.DegreeElevate(p_2 - p_1,1)
end

nu = max(curve_1.nu,curve_2.nu);

[knots_2, knots_1] = KnotToInsert(curve_1.U,curve_2.U);

if knots_2(1) ~= 0
    curve_2.KnotRefine(knots_2,1);
end

if knots_1(1) ~= 0
    curve_1.KnotRefine(knots_1,1);
end

    
U = curve_1.U;

B_1 = curve_1.get_point_cell;
B_2 = curve_2.get_point_cell;

S = cell(curve_1.nu+1,2);
for i = 1:curve_1.nu+1
    S{i,1} = B_1{i};
    S{i,2} = B_2{i};
end

ruled = Geometry('surf',p,curve_1.U,1,[0,0,1,1],S);

end
    
    
    
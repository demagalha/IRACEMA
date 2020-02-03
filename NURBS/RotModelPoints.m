function P = RotModelPoints(Model,face1,face2)

B = Model.get_point_cell;
S = ExtractFacesIndex(Model);



F1 = S{face1}.P;
F2 = S{face2}.P;
clear S;

[sz1,sz2] = size(F1);

P = zeros(sz1*sz2,2);

count = 1;
for i = 1:sz1*sz2
	for j = 1:sz1*sz2
		if isequal([B{F1(i)}(1),B{F1(i)}(2),B{F1(i)}(3)],[B{F2(j)}(1),B{F2(j)}(2),B{F2(j)}(3)])
			P(count,:) = [F1(i),F2(j)];
			count = count + 1;
			break
		end
	end
end

end

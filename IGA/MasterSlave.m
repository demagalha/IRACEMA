function MS = MasterSlave(Model,Model2,face1,face2)

P = Model.get_point_cell; %patch 1;
P2 = Model2.get_point_cell; %patch 2;

SurfCptsIndex = ExtractFacesIndex(Model);
SurfCptsIndex2 = ExtractFacesIndex(Model2);

%suponha que a face 3 do patch 1 se "junte" a face 4 do patch 2 (ie, as 2 faces (faces no espaço paramétrico, superfícies no espaço euclidiano) tem os mesmos pontos de controle)

%%devo achar os "global basis functions number"

%%lembrar que na ID array (ID(i,A)) , i é o grau de liberdade e A o global basis function number
MS = zeros(numel(SurfCptsIndex{face1}.P),2);
count = 1;
for i =1:numel(SurfCptsIndex{face1}.P)
	for j=1:numel(SurfCptsIndex2{face2}.P)
        PP = [P{SurfCptsIndex{face1}.P(i)}(1), P{SurfCptsIndex{face1}.P(i)}(2), P{SurfCptsIndex{face1}.P(i)}(3)];
        PP2 = [P2{SurfCptsIndex{face2}.P(j)}(1), P2{SurfCptsIndex{face2}.P(j)}(2), P2{SurfCptsIndex{face2}.P(j)}(3)];
        
		if isequal(PP,PP2)
			MS(count,:) = [SurfCptsIndex{face1}.P(i),SurfCptsIndex2{face2}.P(j)];
            count = count + 1;
            break
		end
	end
end

end
				
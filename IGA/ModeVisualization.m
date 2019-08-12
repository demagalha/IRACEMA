clear all
close all
clc
load('circplate.mat')
[INN, IEN, ~, ~] = Model.get_connectivity;
 B = Model.get_point_cell;
 x = reshape(B,numel(B),1);
 u = cell(size(x));
 comb = u;
 [sz1 sz2] = size(autovector);
 scaling = 1;
Modos = cell(sz2,1);
for n_modo =1:sz2
    for i=1:size(ID,2)
        u{i} = [autovector(ID(:,i),n_modo); 0]'; 
        comb{i} = x{i} +scaling*u{i};
    end
    tmp = reshape(comb,size(B));
    Modos{n_modo} = Geometry('volume',Model.pu,Model.U,Model.pv,Model.V,Model.pw,Model.W,tmp);
end
% close all
% Modos{4}.plot_geo
% title('Viga Biapoiada, Modo 4','FontSize',23,'FontWeight','bold')
% cd D:\Whirlpool\Simulacao\IGA\Vigas
% cd D:\IRACEMA\NURBS\
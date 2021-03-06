clear all
% close all
clc
load('circplate2.mat')

 B = Model.get_point_cell;
 u = cell(size(B));
 comb = u;
 [sz1 sz2] = size(autovector);
Modos = cell(sz2,1);
 for n_modo = 1:sz2
    for i =1:size(ID,2)
        u{i} = [autovector(3*(i-1) +1,n_modo),autovector(3*(i-1) +2,n_modo),autovector(3*(i-1) +3,n_modo), 0];
        comb{i} = B{i} + u{i};
    end
    Modos{n_modo} = Geometry('volume',Model.pu,Model.U,Model.pv,Model.V,Model.pw,Model.W,comb);
 end
%  subplot(1,2,2)
 Modos{1}.plot_geo
 shading interp
 zlim([0 1])
%  subplot(3,2,4)
%  shading interp
%  zlim([0 1])
%  Modos{2}.plot_geo
%  shading interp
%  zlim([0 1])
%  subplot(3,2,6)
%  Modos{3}.plot_geo
%  zlim([0 1])
%  shading interp
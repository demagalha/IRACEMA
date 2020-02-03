clear all
close all
clc
load('rodk2.mat')
ck2 = coordinate;
fk2 = fun_plot;
load('rodk3.mat')
ck3 = coordinate;
fk3 = fun_plot;
load('rodk4.mat')
ck4 = coordinate;
fk4 = fun_plot;
load('rodp2.mat')
cp2 = coordinate;
fp2 = fun_plot;
load('rodp3.mat')
cp3 = coordinate;
fp3 = fun_plot;
load('rodp4.mat')
cp4 = coordinate;
fp4 = fun_plot;
clear coordinate fun_plot

figure(1)
plot(ck2, fk2,'b', 'LineWidth',2)
hold on
plot(cp2, fp2,'b--','LineWidth',2)
plot(ck3, fk3,'r', 'LineWidth',2)
plot(cp3, fp3,'r--','LineWidth',2)
plot(ck4, fk4,'g', 'LineWidth',2)
plot(cp4, fp4,'g--','LineWidth',2)
legend('IGA, p=2','FEM, p=2','IGA, p=3','FEM, p=3','IGA, p=4','FEM, p=4')
xlabel('n/N','FontWeight','bold','FontSize',23)
ylabel('\omega_n / \omega_{nt}','FontWeight','bold','FontSize',23)
a = annotation('textbox', [0.75, 0.3, 0.1, 0.1], 'String', "NURBS");
a.FontSize = 20;
a.FontWeight = 'bold';
a.EdgeColor = 'white';
b = annotation('textbox', [0.65, 0.65, 0.1 0.1], 'String', "FEA");
b.FontSize = 20;
b.FontWeight = 'bold';
b.EdgeColor = 'white';
set(gca,'FontSize',20)
% figure(2)

% figure(3)
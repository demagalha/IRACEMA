clear all
close all
clc
load ConvergenciaPlacaP5
load ConvergenciaPlacaP4
load ConvergenciaPlacaP3
load ConvergenciaPlacaP2
Placa_p5 = walue;
clear walue

% figure(1)
% semilogx(Placa_p2(7,:),Placa_p2(6,:),'LineWidth',3)
% hold all
% grid on
% semilogx(Placa_p3(7,:),Placa_p3(6,:),'LineWidth',3)
% semilogx(Placa_p4(7,:),Placa_p4(6,:),'LineWidth',3)
% semilogx(Placa_p5(7,:),Placa_p5(6,:),'LineWidth',3)
% legend('p=2','p=3','p=4','p=5')
% xlabel('Number of Elements','FontWeight','bold','FontSize',23)
% ylabel('Time [s]','FontWeight','bold','FontSize',23)
% title('Model Simulation time vs Element Polynomial Order','FontWeight','bold','FontSize',23)
% set(gca,'FontSize',20)

anal = [27; 67.58; 81.57];
ConvP2 = (Placa_p2([1;3;5],:)-anal)./anal;
ConvP3 = (Placa_p3([1;3;5],:)-anal)./anal;
ConvP4 = (Placa_p4([1;3;5],:)-anal)./anal;
ConvP5 = (Placa_p5([1;3;5],:)-anal)./anal;
figure(2)
loglog(Placa_p2(7,:),ConvP2(2,:),'LineWidth',3)
hold all
loglog(Placa_p3(7,:),ConvP3(2,:),'LineWidth',3)
loglog(Placa_p4(7,:),ConvP4(2,:),'LineWidth',3)
loglog(Placa_p5(7,:),ConvP5(2,:),'LineWidth',3)
grid on
legend('p=2','p=3','p=4','p=5')
xlabel('Number of Elements','FontWeight','bold','FontSize',23)
ylabel('Frequency Factor Relative Error [%]','FontWeight','bold','FontSize',23)
title('Mode 3 Frequency Factor vs Element Polynomial Order','FontWeight','bold','FontSize',23)
set(gca,'FontSize',20)

% figure(3)
% loglog(Placa_p2(6,:),ConvP2(3,:),'LineWidth',3)
% hold all
% grid on
% loglog(Placa_p3(6,:),ConvP3(3,:),'LineWidth',3)
% loglog(Placa_p4(6,:),ConvP4(3,:),'LineWidth',3)
% loglog(Placa_p5(6,:),ConvP5(3,:),'LineWidth',3)
% legend('p=2','p=3','p=4','p=5')
% xlabel('Time [s]','FontWeight','bold','FontSize',23)
% ylabel('Frequency Factor Relative Error [%]','FontWeight','bold','FontSize',23)
% title('Mode 5 Frequency Factor vs Time','FontWeight','bold','FontSize',23)
% set(gca,'FontSize',20)

% figure(3)
% loglog(ConvP2(2,:),Placa_p2(6,:),'LineWidth',3)
% hold all
% grid on
% loglog(ConvP3(2,:),Placa_p3(6,:),'LineWidth',3)
% loglog(ConvP4(2,:),Placa_p4(6,:),'LineWidth',3)
% loglog(ConvP5(2,:),Placa_p5(6,:),'LineWidth',3)

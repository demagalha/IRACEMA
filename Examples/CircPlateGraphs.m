clear all
close all
clc
load circplateP2.mat
p2 = walue;
load circplateP3.mat
p3 = walue;
load circplateP4.mat
p4 = walue;
load circplateP5.mat
p5 = walue;
clearvars -except p2 p3 p4 p5

leissa = [10.2158; 21.26; 34.88; 39.771; 51.04; 60.82;];
idx = [1; 2; 4; 6; 7; 9;];

p2_modes = p2(idx,:);
p2_conv = (p2_modes-leissa)./leissa;
p2_conv = [p2_conv ; p2(11,:)];


p3_modes = p3(idx,:);
p3_conv = (p3_modes-leissa)./leissa;
p3_conv = [p3_conv ; p3(11,:)];

p4_modes = p4(idx,:);
p4_conv = (p4_modes-leissa)./leissa;
p4_conv = [p4_conv ; p4(11,:)];

p5_modes = p5(idx,:);
p5_conv = (p5_modes-leissa)./leissa;
p5_conv = [p5_conv ; p5(11,:)];

%% Figures

figure(1)
subplot(2,2,1)
semilogx(p2(12,:),p2_conv(1,:),'LineWidth',3)
hold all
grid on
semilogx(p3(12,:),p3_conv(1,:),'LineWidth',3)
semilogx(p4(12,:),p4_conv(1,:),'LineWidth',3)
semilogx(p5(12,:),p5_conv(1,:),'LineWidth',3)
legend('p=2','p=3','p=4','p=5')
title('Frequency Factor vs Element Polynomial Order','FontWeight','bold','FontSize',23)
ylabel('Frequency Factor Relative Error ','FontWeight','bold','FontSize',23)
xlabel('Number of Elements','FontWeight','bold','FontSize',23)
set(gca,'FontSize',20,'FontWeight','bold')

% figure(2)
% subplot(3,1,2)
% semilogx(p2(12,:),p2_conv(2,:),'LineWidth',3)
% hold all
% grid on
% semilogx(p3(12,:),p3_conv(2,:),'LineWidth',3)
% semilogx(p4(12,:),p4_conv(2,:),'LineWidth',3)
% semilogx(p5(12,:),p5_conv(2,:),'LineWidth',3)
% legend('p=2','p=3','p=4','p=5')
% title('Frequency Factor vs Element Polynomial Order','FontWeight','bold','FontSize',23)
% ylabel('Frequency Factor Relative Error ','FontWeight','bold','FontSize',23)
% xlabel('Number of Elements','FontWeight','bold','FontSize',23)
% set(gca,'FontSize',20,'FontWeight','bold')

% figure(3)
% subplot(3,1,3)
% semilogx(p2(12,:),p2_conv(3,:),'LineWidth',3)
% hold all
% grid on
% semilogx(p3(12,:),p3_conv(3,:),'LineWidth',3)
% semilogx(p4(12,:),p4_conv(3,:),'LineWidth',3)
% semilogx(p5(12,:),p5_conv(3,:),'LineWidth',3)
% legend('p=2','p=3','p=4','p=5')
% title('Frequency Factor vs Element Polynomial Order','FontWeight','bold','FontSize',23)
% ylabel('Frequency Factor Relative Error ','FontWeight','bold','FontSize',23)
% xlabel('Number of Elements','FontWeight','bold','FontSize',23)
% set(gca,'FontSize',20,'FontWeight','bold')

% figure(4)
% semilogx(p2(12,:),p2_conv(4,:),'LineWidth',3)
% hold all
% grid on
% semilogx(p3(12,:),p3_conv(4,:),'LineWidth',3)
% semilogx(p4(12,:),p4_conv(4,:),'LineWidth',3)
% semilogx(p5(12,:),p5_conv(4,:),'LineWidth',3)
% legend('p=2','p=3','p=4','p=5')
% title('Frequency Factor vs Element Polynomial Order','FontWeight','bold','FontSize',23)
% ylabel('Frequency Factor Relative Error ','FontWeight','bold','FontSize',23)
% xlabel('Number of Elements','FontWeight','bold','FontSize',23)
% set(gca,'FontSize',20,'FontWeight','bold')

% figure(5)
% semilogx(p2(12,:),p2_conv(5,:),'LineWidth',3)
% hold all
% grid on
% semilogx(p3(12,:),p3_conv(5,:),'LineWidth',3)
% semilogx(p4(12,:),p4_conv(5,:),'LineWidth',3)
% semilogx(p5(12,:),p5_conv(5,:),'LineWidth',3)
% legend('p=2','p=3','p=4','p=5')
% title('Frequency Factor vs Element Polynomial Order','FontWeight','bold','FontSize',23)
% ylabel('Frequency Factor Relative Error ','FontWeight','bold','FontSize',23)
% xlabel('Number of Elements','FontWeight','bold','FontSize',23)
% set(gca,'FontSize',20,'FontWeight','bold')

% figure(6)
% semilogx(p2(12,:),p2_conv(6,:),'LineWidth',3)
% hold all
% grid on
% semilogx(p3(12,:),p3_conv(6,:),'LineWidth',3)
% semilogx(p4(12,:),p4_conv(6,:),'LineWidth',3)
% semilogx(p5(12,:),p5_conv(6,:),'LineWidth',3)
% legend('p=2','p=3','p=4','p=5')
% title('Frequency Factor vs Element Polynomial Order','FontWeight','bold','FontSize',23)
% ylabel('Frequency Factor Relative Error ','FontWeight','bold','FontSize',23)
% xlabel('Number of Elements','FontWeight','bold','FontSize',23)
% set(gca,'FontSize',20,'FontWeight','bold')

% figure(7)
% semilogx(p2(11,:),p2_conv(1,:),'LineWidth',3)
% hold all
% grid on
% semilogx(p3(11,:),p3_conv(1,:),'LineWidth',3)
% semilogx(p4(11,:),p4_conv(1,:),'LineWidth',3)
% semilogx(p5(11,:),p5_conv(1,:),'LineWidth',3)
% legend('p=2','p=3','p=4','p=5')
% title('Frequency Factor vs Element Polynomial Order','FontWeight','bold','FontSize',23)
% ylabel('Frequency Factor Relative Error ','FontWeight','bold','FontSize',23)
% xlabel('Number of Elements','FontWeight','bold','FontSize',23)
% set(gca,'FontSize',20,'FontWeight','bold')

% figure(8)
% semilogx(p2(11,:),p2_conv(2,:),'LineWidth',3)
% hold all
% grid on
% semilogx(p3(11,:),p3_conv(2,:),'LineWidth',3)
% semilogx(p4(11,:),p4_conv(2,:),'LineWidth',3)
% semilogx(p5(11,:),p5_conv(2,:),'LineWidth',3)
% legend('p=2','p=3','p=4','p=5')
% title('Frequency Factor vs Element Polynomial Order','FontWeight','bold','FontSize',23)
% ylabel('Frequency Factor Relative Error ','FontWeight','bold','FontSize',23)
% xlabel('Number of Elements','FontWeight','bold','FontSize',23)
% set(gca,'FontSize',20,'FontWeight','bold')

% figure(9)
% semilogx(p2(11,:),p2_conv(3,:),'LineWidth',3)
% hold all
% grid on
% semilogx(p3(11,:),p3_conv(3,:),'LineWidth',3)
% semilogx(p4(11,:),p4_conv(3,:),'LineWidth',3)
% semilogx(p5(11,:),p5_conv(3,:),'LineWidth',3)
% legend('p=2','p=3','p=4','p=5')
% title('Frequency Factor vs Element Polynomial Order','FontWeight','bold','FontSize',23)
% ylabel('Frequency Factor Relative Error ','FontWeight','bold','FontSize',23)
% xlabel('Number of Elements','FontWeight','bold','FontSize',23)
% set(gca,'FontSize',20,'FontWeight','bold')

% figure(10)
% semilogx(p2(11,:),p2_conv(4,:),'LineWidth',3)
% hold all
% grid on
% semilogx(p3(11,:),p3_conv(4,:),'LineWidth',3)
% semilogx(p4(11,:),p4_conv(4,:),'LineWidth',3)
% semilogx(p5(11,:),p5_conv(4,:),'LineWidth',3)
% legend('p=2','p=3','p=4','p=5')
% title('Frequency Factor vs Element Polynomial Order','FontWeight','bold','FontSize',23)
% ylabel('Frequency Factor Relative Error ','FontWeight','bold','FontSize',23)
% xlabel('Number of Elements','FontWeight','bold','FontSize',23)
% set(gca,'FontSize',20,'FontWeight','bold')

% figure(11)
% semilogx(p2(11,:),p2_conv(5,:),'LineWidth',3)
% hold all
% grid on
% semilogx(p3(11,:),p3_conv(5,:),'LineWidth',3)
% semilogx(p4(11,:),p4_conv(5,:),'LineWidth',3)
% semilogx(p5(11,:),p5_conv(5,:),'LineWidth',3)
% legend('p=2','p=3','p=4','p=5')
% title('Frequency Factor vs Element Polynomial Order','FontWeight','bold','FontSize',23)
% ylabel('Frequency Factor Relative Error ','FontWeight','bold','FontSize',23)
% xlabel('Number of Elements','FontWeight','bold','FontSize',23)
% set(gca,'FontSize',20,'FontWeight','bold')
% 
% figure(12)
% semilogx(p2(11,:),p2_conv(6,:),'LineWidth',3)
% hold all
% grid on
% semilogx(p3(11,:),p3_conv(6,:),'LineWidth',3)
% semilogx(p4(11,:),p4_conv(6,:),'LineWidth',3)
% semilogx(p5(11,:),p5_conv(6,:),'LineWidth',3)
% legend('p=2','p=3','p=4','p=5')
% title('Frequency Factor vs Element Polynomial Order','FontWeight','bold','FontSize',23)
% ylabel('Frequency Factor Relative Error ','FontWeight','bold','FontSize',23)
% xlabel('Number of Elements','FontWeight','bold','FontSize',23)
% set(gca,'FontSize',20,'FontWeight','bold')

% figure(13)
subplot(2,2,3)
semilogx(p2(12,:),p2(11,:),'LineWidth',3)
hold all
grid on
semilogx(p3(12,:),p3(11,:),'LineWidth',3)
semilogx(p4(12,:),p4(11,:),'LineWidth',3)
semilogx(p5(12,:),p5(11,:),'LineWidth',3)
legend('p=2','p=3','p=4','p=5')
title('Simulation Time vs Element Polynomial Order','FontWeight','bold','FontSize',23)
ylabel('Time [s]','FontWeight','bold','FontSize',23)
xlabel('Number of Elements','FontWeight','bold','FontSize',23)
set(gca,'FontSize',20,'FontWeight','bold')
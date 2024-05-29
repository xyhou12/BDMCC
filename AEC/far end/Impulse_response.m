clc
close all;
clear all;

load car_echo_path.mat
load h_shift.mat h_shift

figure
plot(1:256,h,'-','LineWidth',1.5)

grid on

xlabel('Numbers of Taps')
ylabel('Amplitude')

set(get(gca,'XLabel'),'Fontname', 'Times New Roman','FontSize',11);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'Fontname', 'Times New Roman','FontSize',11);
set(get(gca,'TITLE'),'Fontname', 'Times New Roman','FontSize',11);
set(gca,'Fontname', 'Times New Roman','fontsize',11);
% set(l1,'Fontname', 'Times New Roman','FontSize',11)
xlim([0,255])
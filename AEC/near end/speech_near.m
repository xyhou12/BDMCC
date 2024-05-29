clc 
close all
clear all

iter = 160000;
load near_end.mat

hold on
plot(near_end)
grid on
xlabel('Samples')
ylabel('Magnitude')

set(get(gca,'XLabel'),'Fontname', 'Times New Roman','FontSize',11);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'Fontname', 'Times New Roman','FontSize',11);
set(get(gca,'TITLE'),'Fontname', 'Times New Roman','FontSize',11);
set(gca,'Fontname', 'Times New Roman','fontsize',11);
% l1 = legend('MCC(\sigma = 1, \mu = 0.005)','NMCC(\sigma = 1, \mu = 0.25)','APMCC({\itP} = 2, \sigma = 1, \mu = 0.15)',['BDMCC({\itM} = 2, \sigma_{1} = \sigma_{2} = 4, ' newline '\lambda_{\itq} = 0.999996, \lambda_{\itz} = 0.9999, \lambda_{\itv} = 0.99999)'],'Far-end speech');
% l1 = legend('MCC','NMCC','APMCC','BDMCC');
% set(l1,'Fontname', 'Times New Roman','FontSize',11)

% plot(1:iter,10*log10(MSD_th_1(1:5000)),'b-')
% plot(1:iter,10*log10(MSD_th_2(1:5000)),'b-')
% plot(1:iter,10*log10(MSD_th_5(1:5000)),'b-')
% % plot(1:iter,10*log10(MSD_th_10(1:5000)),'b-')
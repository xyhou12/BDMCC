% NSAF
clc; % 擦去一页命令窗口，光标回屏幕左上角
clear all; % 从工作空间清除所有变量和函数
close all; % 关闭所有窗口
tic
% ………………运行参数……………………
run = 1; % 独立运行次数
N = 400000; % 采样点数/迭代次数
L = 800;%待识别的系统的阶数
N = 4; % 子带数目
% load filter_bank_4.mat % hk：分析滤波器, gk：综合滤波器，子带数4
% hk = Hkn;
% ………………系统初始化……………………
load speech40k.mat  %语音信
load blocksystem.mat W1;
unknown_w = W1;
unknown_w1 = W1/norm(W1);
% unknown_w2 = -unknown_w1;

delta = 1;
delta1 = 1;
a_hat = zeros(M,L+P+1);
sigma_z_n = 8;
sigma_v_n = 1;
sigma_q_n = 0;
sigma_a = 1;
sigma = 1;
D = zeros(1,M);

sigma_v_hat = 0;
lambda_q = 0.999999;
lambda_z = 0.999;
lambda_v = 0.999999;

delta_q = 1e-10;
mu1 = 1.2; % SSAF固定步长
MSD3 = zeros(1,N); % 均方权值偏差
%………………算法……………………
for g = 1:run
    g
    % ………………算法初始化……………………
    w1 = zeros(L,1); % 权值初始化
    % ………………输入信号……………………
    xx = CC';
    U = 30*xx(1:N);
    var_u = var(U);
    delta = 10*var_u;
    delta_mao = 100*ones(sub_N,1);
    alpha = zeros(sub_N,1);
    %………………期望信号……………………
%     d = filter(unknown_w1,1,U); % 未知系统输出信号
    d1 = filter(unknown_w1,1,U(1:N/2)); % 未知系统输出信号
    d2 = filter(unknown_w2,1,U(N/2+1:N));
    d = [d1 d2];
%     d1 = filter(c1,1,U(1:N/3)); % 未知系统输出信号
%     d2 = filter(c2,1,U((N/3+1):2*N/3));
%     d3 = filter(c3,1,U((2*N/3+1):N));
%     d = [d1 d2 d3];
    D = awgn(d,30); % 加入30db高斯白噪声
    var_v = var(D-d);
    var_v_sub = var_v/sub_N;
    %………………子带分解……………………
    sub_U = zeros(sub_N,N); %输入信号
    sub_D = zeros(sub_N,N/sub_N); %子带期望信号
    for j1 = 1:sub_N
        sub_U(j1,:) = filter(hk(j1,:),1,U); % 输入信号分解   % hk：分析滤波器
        sub_d = filter(hk(j1,:),1,D); % 期望信号分解
        sub_D(j1,:) = sub_d(1:sub_N:N); % 抽取
    end
    %………………更新过程………………
    sub_x = zeros(L,sub_N); %子带输入信号
    sub_e1 = zeros(sub_N,1); %子带误差信号
    sub_DD = zeros(sub_N,1);
    for i = L/sub_N+1:N/sub_N-1
        sub_x = sub_U(:,i*sub_N-(sub_N-1):-1:i*sub_N-(sub_N-1)-(L-1))';
        sub_DD = sub_D(:,i);
        %——————NSAF————————
        sub_e1 = sub_DD - sub_x'*w1;
        sum_xe = zeros(L,1);
        for ii = 1:sub_N
            alpha(ii) = ( delta_mao(ii) + delta_q ) / ( ( delta_mao(ii) + delta_q ) * (sub_x(:,ii)'*sub_x(:,ii)) + var_v_sub  );
            delta_mao(ii) = ( 1 - ( alpha(ii) * (sub_x(:,ii)'*sub_x(:,ii)) ) / L  ) * ( delta_mao(ii) + delta_q  );
        end
        for ii1 = 1:sub_N
            sum_xe = sum_xe + alpha(ii1)*sub_e1(ii1)*sub_x(:,ii1);
        end
        w1 = w1 + (1/sub_N) * sum_xe; % 滤波器权系数更新
%         MSD1(i*sub_N) = MSD1(i*sub_N) + (norm(w1 - unknown_w1)/norm(unknown_w1))^2;
        if i*sub_N <= N/2
            MSD3(i*sub_N) = MSD3(i*sub_N) + (norm(w1 - unknown_w1)/norm(unknown_w1))^2;
        else
            MSD3(i*sub_N) = MSD3(i*sub_N) + (norm(w1 - unknown_w2)/norm(unknown_w2))^2;
        end
%         if i*sub_N <= N/3
%             MSD3(i*sub_N) = MSD3(i*sub_N) + (norm(w1 - c1)/norm(c1))^2;
%         elseif i*sub_N >= (N/3+1) && i*sub_N <= (2*N/3)
%             MSD3(i*sub_N) = MSD3(i*sub_N) + (norm(w1 - c2)/norm(c2))^2;
%         else  i*sub_N >= 2*N/3+1
%             MSD3(i*sub_N) = MSD3(i*sub_N) + (norm(w1 - c3)/norm(c3))^2;
%         end
        MSD3(i*sub_N+(1:sub_N-1)) = MSD3(i*sub_N); % 插值
    end
end
%………………结果输出MSD……………………
MSD3 = MSD3/run;
save BSAF.mat MSD3;
figure(2)
plot(1:N,10*log10(MSD3),'b');
xlabel('迭代次数');
ylabel('MSD(dB)');
legend('NSAF');
title('MSD(dB)');
grid on;
toc
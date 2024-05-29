% NSAF
clc; % ��ȥһҳ����ڣ�������Ļ���Ͻ�
clear all; % �ӹ����ռ�������б����ͺ���
close all; % �ر����д���
tic
% ���������������в�������������������
run = 1; % �������д���
N = 400000; % ��������/��������
L = 800;%��ʶ���ϵͳ�Ľ���
N = 4; % �Ӵ���Ŀ
% load filter_bank_4.mat % hk�������˲���, gk���ۺ��˲������Ӵ���4
% hk = Hkn;
% ������������ϵͳ��ʼ������������������
load speech40k.mat  %������
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
mu1 = 1.2; % SSAF�̶�����
MSD3 = zeros(1,N); % ����Ȩֵƫ��
%�������������㷨����������������
for g = 1:run
    g
    % �������������㷨��ʼ������������������
    w1 = zeros(L,1); % Ȩֵ��ʼ��
    % �����������������źš���������������
    xx = CC';
    U = 30*xx(1:N);
    var_u = var(U);
    delta = 10*var_u;
    delta_mao = 100*ones(sub_N,1);
    alpha = zeros(sub_N,1);
    %�����������������źš���������������
%     d = filter(unknown_w1,1,U); % δ֪ϵͳ����ź�
    d1 = filter(unknown_w1,1,U(1:N/2)); % δ֪ϵͳ����ź�
    d2 = filter(unknown_w2,1,U(N/2+1:N));
    d = [d1 d2];
%     d1 = filter(c1,1,U(1:N/3)); % δ֪ϵͳ����ź�
%     d2 = filter(c2,1,U((N/3+1):2*N/3));
%     d3 = filter(c3,1,U((2*N/3+1):N));
%     d = [d1 d2 d3];
    D = awgn(d,30); % ����30db��˹������
    var_v = var(D-d);
    var_v_sub = var_v/sub_N;
    %�������������Ӵ��ֽ⡭��������������
    sub_U = zeros(sub_N,N); %�����ź�
    sub_D = zeros(sub_N,N/sub_N); %�Ӵ������ź�
    for j1 = 1:sub_N
        sub_U(j1,:) = filter(hk(j1,:),1,U); % �����źŷֽ�   % hk�������˲���
        sub_d = filter(hk(j1,:),1,D); % �����źŷֽ�
        sub_D(j1,:) = sub_d(1:sub_N:N); % ��ȡ
    end
    %���������������¹��̡�����������
    sub_x = zeros(L,sub_N); %�Ӵ������ź�
    sub_e1 = zeros(sub_N,1); %�Ӵ�����ź�
    sub_DD = zeros(sub_N,1);
    for i = L/sub_N+1:N/sub_N-1
        sub_x = sub_U(:,i*sub_N-(sub_N-1):-1:i*sub_N-(sub_N-1)-(L-1))';
        sub_DD = sub_D(:,i);
        %������������NSAF����������������
        sub_e1 = sub_DD - sub_x'*w1;
        sum_xe = zeros(L,1);
        for ii = 1:sub_N
            alpha(ii) = ( delta_mao(ii) + delta_q ) / ( ( delta_mao(ii) + delta_q ) * (sub_x(:,ii)'*sub_x(:,ii)) + var_v_sub  );
            delta_mao(ii) = ( 1 - ( alpha(ii) * (sub_x(:,ii)'*sub_x(:,ii)) ) / L  ) * ( delta_mao(ii) + delta_q  );
        end
        for ii1 = 1:sub_N
            sum_xe = sum_xe + alpha(ii1)*sub_e1(ii1)*sub_x(:,ii1);
        end
        w1 = w1 + (1/sub_N) * sum_xe; % �˲���Ȩϵ������
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
        MSD3(i*sub_N+(1:sub_N-1)) = MSD3(i*sub_N); % ��ֵ
    end
end
%������������������MSD����������������
MSD3 = MSD3/run;
save BSAF.mat MSD3;
figure(2)
plot(1:N,10*log10(MSD3),'b');
xlabel('��������');
ylabel('MSD(dB)');
legend('NSAF');
title('MSD(dB)');
grid on;
toc
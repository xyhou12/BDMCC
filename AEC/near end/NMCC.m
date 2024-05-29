% clc;
% clear;
% close all;

L = 2048;
mu = 2.2;
P = 2;   %投影阶数
sigma = 0.3;

load speech.mat
x = 5*xx(1:160000)';
near_end = [zeros(1,40000) 5*xx(40001:60000)' zeros(1,60000) 5*xx(40001:60000)' zeros(1,20000)];%%近端语音

% plot(near_end)

iter = length(x);
deta = 20*mean(x.^2);
% deta = 0.01;

%%回波路径
% load car_echo_path.mat
% load h_shift.mat h_shift

load channe.mat
h = gL;
h_shift = h;
% h_shift = h;
% h = h_shift;
d1 = filter(h,1,x(1:iter/2));
d2 = zeros(1,iter/2);
for ii1 = 1:iter/2
    d2(ii1) = x(iter/2+ii1:-1:iter/2+ii1-L+1)*h_shift;
end
d = [d1 d2];
var_noise = 1e-3;
noise = sqrt(var_noise)*randn(1,iter);
sas = stblrnd(1.2,0,1/30,0,1,iter);
% sas = rasd(iter, 1.0, 0, 0.1, 0);
% plot(sas)
d_noise = d + noise + sas + near_end;

w = zeros(L,1);
E = zeros(P,1);
DD = zeros(P,1);
r = zeros(1,iter);
X = zeros(L,P);
NMSD = zeros(1,iter);
for k = L+P:iter
    x_input = x(k:-1:k-L+1)';
    y = x_input'*w;
    e = d_noise(k) - y;
    Xs = x_input*exp(-e^2/(2*sigma^2))*e/(x_input'*x_input + deta);
    w = w + mu*Xs;
    if k <= iter/2
        NMSD(k) = NMSD(k) + (norm(w - h)/norm(h))^2;
    else
        NMSD(k) = NMSD(k) + (norm(w - h_shift)/norm(h_shift))^2;
    end
end

%%NMSD
% figure
plot(x)
hold on
plot(10*log10(NMSD),'b-','linewidth',1.5)
grid on
xlim([1 iter])
set(gca,'FontSize',12);

% 
% NMSD_NLMS = NMSD;
% save NLMS.mat NMSD_NLMS
NMSD_NMCC = NMSD;
save NMCC.mat NMSD_NMCC
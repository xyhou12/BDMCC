% clc;
% clear;
% close all;


L = 2048;
mu = 0.4;
P = 4;   %投影阶数
sigma = 0.3;

load speech.mat
x = 5*xx(1:160000)';
% [x,PS] = mapstd(x);
near_end = [zeros(1,40000) 5*xx(40001:60000)' zeros(1,60000) 5*xx(40001:60000)' zeros(1,20000)];%%近端语音

% plot(near_end)

iter = length(x);
deta = 20*mean(x.^2);
% deta = 1;

%%回波路径
% load car_echo_path.mat
% load h_shift.mat h_shift
load channe.mat
h = gL;
h_shift = h;

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
d_noise = d + noise + sas;

w = zeros(L,1);
E = zeros(P,1);
DD = zeros(P,1);
r = zeros(1,iter);
X = zeros(L,P);
NMSD = zeros(1,iter);
for k = L+P:iter
    DD = d_noise(k:-1:k-P+1)';
    for kk1 = 1:P
        X(:, kk1) = x(k-kk1+1:-1:k-kk1+1-L+1)';
    end
    Y = X'*w;
    y = Y(1);
    E = DD - Y;
    r(k) = d(k) - y;
%     Xs = X*sign(E);
%     Xs = X*diag(exp(-E.^2/(2*sigma^2)))*E;
    w = w + mu*X*inv(X'*X +deta*eye(P))*diag(exp(-E.^2/(2*sigma^2)))*E;
    if k <= iter/2
        NMSD(k) = NMSD(k) + (norm(w - h)/norm(h))^2;
    else
        NMSD(k) = NMSD(k) + (norm(w - h_shift)/norm(h_shift))^2;
    end
end

%%NMSD
% figure
% plot(x)
hold on
plot(10*log10(NMSD),'b-','linewidth',1.5)
grid on
xlim([1 iter])
set(gca,'FontSize',12);

% NMSD_APLMS = NMSD;
% save APLMS.mat NMSD_APLMS 
NMSD_APMCC = NMSD;
save APMCC.mat NMSD_APMCC 
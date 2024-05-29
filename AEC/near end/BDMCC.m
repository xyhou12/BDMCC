% clc
% clear
% close all

L = 2048;
M = 2;

load speech.mat
x = 5*xx(1:160000)';
near_end = [zeros(1,40000) 5*xx(40001:60000)' zeros(1,60000) 5*xx(40001:60000)' zeros(1,20000)];%%½ü¶ËÓïÒô

iter = length(x);

%%»Ø²¨Â·¾¶
% load car_echo_path.mat
% load h_shift.mat h_shift
load channe.mat
h = gL;
h_shift = h;
% h = h_shift;
% h = h(1:64);

d1 = filter(h,1,x(1:iter/2));
d2 = zeros(1,iter/2);
for ii1 = 1:iter/2
    d2(ii1) = x(iter/2+ii1:-1:iter/2+ii1-L+1)*h_shift;
end
d = [d1 d2];

var_noise = 1e-3;
noise = sqrt(var_noise)*randn(1,iter);
sas = stblrnd(1.2,0,1/30,0,1,iter);
d_noise = d + noise + sas + near_end;

w = zeros(L,L+M+1);
r = zeros(1,iter);
X = zeros(L,M);
NMSD = zeros(1,iter);

delta = 0.3;
delta1 = 1;
a_hat = zeros(M,L+M+1);
sigma_z_n = 1;
sigma_v_n = 1;
sigma_q_n = 0;

sigma_w = 0;
sigma_a = 1;
sigma = 1;
D = zeros(1,M);

sigma_v_hat = 0;
lambda_q = 0.99999;
lambda_z = 0.99999;
lambda_v = 0.999999;
for n = L+M+1:iter       
    
    x_a = x(n-1:-1:n-M)';
    x_input = x(n:-1:n-L+1)';
    X = [x(n-1:-1:n-L)',X(:,1:end-1)];
        
    y = x_input'*w(:,n);
   
    e_a = a_hat(:,n) - a_hat(:,n-1);
    e_x = x(n) - x_a'*a_hat(:,n);
    K = exp(-e_x^2/sigma_z_n/2/delta1^2)/exp(-norm(e_a,2)^2/sigma_a/2/delta1^2);
    alpha = K*sigma_a/(sigma_a*norm(x_a,2)^2 + sigma_z_n);
    sigma_a = (1 - alpha*norm(x_a,2)^2/M)*sigma_a;
    
    a_hat(:,n+1) = a_hat(:,n) + alpha*e_x*x_a;
    
    x_hat = x_input - X*a_hat(:,n+1);
%     d_hat = y - D*a_hat(:,n+1);
    
    d_hat = d_noise(n) - D*a_hat(:,n+1);
    
    beta = sigma + sigma_q_n;
       
    e_hat = d_hat - x_hat'*w(:,n);
%     sigma_v_n = (1 + a_hat(:,n+1)'*a_hat(:,n+1))*sigma_v;
%     sigma_v_n = q*sigma_v;
    sigma_v_hat = sigma_v_hat + sigma_v_n;
    e_w =  w(:,n) - w(:,n-1);
    K = exp(-e_hat^2/sigma_v_n/2/delta^2)/exp(-norm(e_w,2)^2/2/beta/delta^2);
    alpha = K*beta/(beta*norm(x_hat,2)^2 + sigma_v_n);
    sigma = (1 - alpha*norm(x_hat,2)^2/L)*beta;
    
    w(:,n+1) = w(:,n) + alpha*e_hat*x_hat;
    D = [d_noise(n),D(1:end-1)];
    
    sigma_w = lambda_q*sigma_w + (1 - lambda_q)*norm(w(:,n+1) - w(:,n),2)^2;
    sigma_q_n = sigma_w/L; 
    sigma_z_n = lambda_z*sigma_z_n  + (1 - lambda_z)*(x_a'*(a_hat(:,n+1) - a_hat(:,n)))^2;
    sigma_v_n = lambda_v*sigma_v_n  + (1 - lambda_v)*(x_hat'*(w(:,n+1) - w(:,n)))^2;
%     NMSD(n) = NMSD(n) + (norm(w(:,n+1) - h)/norm(h))^2;
    if n <= iter/2
        NMSD(n) = NMSD(n) + (norm(w(:,n+1) - h)/norm(h))^2;
    else
        NMSD(n) = NMSD(n) + (norm(w(:,n+1) - h_shift)/norm(h_shift))^2;
    end
end

%%NMSD
% figure
% plot(x)
hold on
plot(10*log10(NMSD),'-','linewidth',1.5)
grid on
xlim([1 iter])
set(gca,'FontSize',12);

%%ERLE

NMSD_BDMCC = NMSD;
save BDMCC.mat NMSD_BDMCC
% NMSD_BDLMS = NMSD;
% save BDLMS.mat NMSD_BDLMS
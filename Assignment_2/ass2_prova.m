%% Comparison of the exponential and running mean for random walk model

close all; clear all; clc;

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 1.1) Generate a true trajectory 𝑋𝑖 using the random walk model 

var_w = 12;
var_eta = 9;

n = 300;
% n = 3000;
x(1) = 10;

% Random Numbers from Normal Distribution with Zero Mean and Variance
a_1 = sqrt(var_w);
a_2 = sqrt(var_eta);

w = a_1.*randn(n,1);
eta = a_2.*randn(n,1);

for i = 2:n
    x(i) = x(i-1) + w(i); % generated trajectory RWM
end

for i=1:n
    z(i) = x(i) + eta(i); % Generate measurements 𝑧𝑖 of the process 𝑋𝑖 
end


%% 2) Identify  𝜎𝑤2  and  𝜎𝜂2  using  identification  method

E_v_sq_sum = [];
for i = 2:n
    E_v_sq_sum(i-1) = ( w(i) + eta(i) - eta(i-1) )^2;
end

E_v_sq = 1/(n-1) *sum(E_v_sq_sum);

E_rho_sq_sum = [];
for i=3:n
     E_rho_sq_sum(i-2) = ( w(i) + w(i-1) + eta(i) - eta(i-2) )^2;
end

E_rho_sq = 1/(n-2) *sum(E_rho_sq_sum);

% try = inv([1 2; 2 2])*[E_v_sq; E_rho_sq]; to delete

syms A B

% A = sigma_w^2
% B = sigma_eta^2

eqns = [ A -  E_v_sq + 2*B == 0,...
          2*B + 2*A - E_rho_sq == 0  ];
vars = [A B];

[a, b] = solve(eqns,vars);
sigma_w2 = double(a)
sigma_eta2 = double(b)

%% 3) Determine optimal smoothing coefficient in exponential smoothing 

csi = sigma_w2/sigma_eta2; 
alpha = (-csi + sqrt(csi^2 + 4*csi))/2; % correct bc should be between 0,1

%% 4) Perform exponential smoothing with the determined smoothing coefficient

x_hat(1) = 10;

for i = 2:n
    x_hat(i) = alpha*z(i) + (1-alpha)*x_hat(i-1);
end

figure(1)
plot(x, 'c', 'LineWidth', 1.2)
hold on
plot(x_hat, 'k', 'LineWidth', 1.2)
plot(z, 'm', 'LineWidth', 1.2)
grid on; grid minor
xlabel('Steps', 'FontSize', 30)
ylabel('Data', 'FontSize', 30)
legend('True Data', 'Smoothed Data', 'Measurements', 'FontSize', 30) 
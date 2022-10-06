%% Determining and removing drawbacks of exponential and running mean

close all; clear; clc;

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Backward exponential smoothing

n_3 = 300; % size of trajectory
incond = 10; % initial condition
x_n(1) = incond;

sigma_w_n = 28^2; % variance noise
sigma_eta_n = 97^2; % variance of noise measurement

a_1_n = sqrt(sigma_w_n);
a_2_n = sqrt(sigma_eta_n);

w_n = a_1_n.*randn(n_3,1);
eta_n = a_2_n.*randn(n_3,1);

for i = 2:n_3
    x_n(i) = x_n(i-1) + w_n(i); % generated trajectory RWM
end

for i=1:n_3
    z_n(i) = x_n(i) + eta_n(i); % Generate measurements 𝑧𝑖 of the process 𝑋𝑖 
end

csi_n = sigma_w_n / sigma_eta_n;
alpha_n = (-csi_n + sqrt(csi_n^2 + 4*csi_n))/2; % correct bc should be between 0,1

% Window size M

M = round((2-alpha_n)/alpha_n); % 7

% Running mean (last measurements are used)
j = (M-1)/2;

x_hat_run = zeros(n_3,1);

x_hat_run(1:j,1) = sum(z_n(1:j))/3;
x_hat_run((n_3-j+1):n_3,1) = sum(z_n((n_3-j+1):n_3))/3;

for i = j+1:n_3-j-1
    x_hat_run(i) = 1/M * (z_n(i-3)+ z_n(i-2) + z_n(i-1) + z_n(i) + ...
    z_n(i+1) + z_n(i+2) + z_n(i+3));
end


% Forward Exponential Estimates
x_hat_forw(1) = incond;

for i = 2:n_3
    x_hat_forw(i) = x_hat_forw(i - 1) + alpha_n*(z_n(i) - x_hat_forw(i - 1));
end

% Apply Backward Smoothing
x_hat_back(n_3) = x_hat_forw(end);

for i = (n_3 - 1):-1:1
    x_hat_back(i) = x_hat_back(i + 1) + alpha_n*(x_hat_forw(i) - x_hat_back(i + 1));
end

figure(1)
hold on
plot(x_n, 'k', 'LineWidth', 1.2)
plot(z_n, 'y', 'LineWidth', 1.2)
plot(x_hat_run, 'b', 'LineWidth', 1.2)
plot(x_hat_back, 'r', 'LineWidth', 1.2)
grid on; grid minor
xlabel('Steps', 'FontSize', 30)
ylabel('Data', 'FontSize', 30)
legend('Trajectory', 'Measuraments', 'Running Mean', 'Backward Exponential Smoothing', 'FontSize', 30)

%% Calculate the indicators

% Deviation Indicators

dev_ind_back = [];
dev_ind_run = [];

for i = 1:n_3

    dev_ind_back(i) = (z_n(i) - x_hat_back(i))^2;
    dev_ind_run(i) = (z_n(i) - x_hat_run(i))^2;

end

dev_ind_back = sum(dev_ind_back)
dev_ind_run = sum(dev_ind_run)

% Variablility Indicators

var_ind_back = [];
var_ind_run = [];

for i = 1:(n_3 - 2)

    var_ind_back(i) = (x_hat_back(i + 2) - 2*x_hat_back(i + 1) + x_hat_back(i))^2;
    var_ind_run(i) = (x_hat_run(i + 2) - 2*x_hat_run(i + 1) + x_hat_run(i))^2;

end

var_ind_back = sum(var_ind_back)
var_ind_run = sum(var_ind_run)

%% Part 2. Drawbacks of running mean

%% 1. Generate a true trajectory 𝑋𝑖 of an object motion disturbed by normally...
% distributed random acceleration 

% n = n_3 = 300

x(1) = 5;
v(1) = 0;
t = 0.1;

sigma_a2 = 10; % variance of noise
sigma_eta2 = 500;

a = sqrt(sigma_a2).*randn(n_3,1);
eta = sqrt(sigma_eta2).*randn(n_3,1);

for i = 2:n_3
    x(i) = x(i - 1) + v(i - 1)*t + (a(i - 1)*t^2)/2;
    v(i) = v(i - 1) + a(i - 1)*t;
end

z = [];

for i = 1:n_3
    z(i) = x(i) + eta(i);
end

M = 80;
x_hat_run_new = movmean(z, M);

figure(2)
plot(x, 'r', 'LineWidth', 1.2)
hold on
plot(z, 'k', 'LineWidth', 1.2)
plot(x_hat_run_new, 'b', 'LineWidth', 1.2)
grid on; grid minor
xlabel('Steps', 'FontSize', 30)
ylabel('Data', 'FontSize', 30)
legend('Trajectory', 'Measuraments', 'Running Mean', 'FontSize', 30)


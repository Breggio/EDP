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

% Running mean
M_guess = 80;
x_hat_run_new = movmean(z, M_guess);

% Forward mean
alpha_guess = 0.005;
x_hat_forw_n(1) = incond;

for i = 2:n_3
    x_hat_forw_n(i) = x_hat_forw_n(i - 1) + alpha_guess*(z_n(i) - x_hat_forw_n(i - 1));
end

figure(2)
plot(x, 'r', 'LineWidth', 1.2)
hold on
plot(z, 'k', 'LineWidth', 1.2)
plot(x_hat_run_new, 'b', 'LineWidth', 1.2)
grid on; grid minor
xlabel('Steps', 'FontSize', 30)
ylabel('Data', 'FontSize', 30)
legend('Trajectory', 'Measurements', 'Running Mean', 'FontSize', 30)

figure(3)
plot(x, 'r', 'LineWidth', 1.2)
hold on
plot(z, 'k', 'LineWidth', 1.2)
plot(x_hat_forw_n, 'b', 'LineWidth', 1.2)
grid on; grid minor
xlabel('Steps', 'FontSize', 30)
ylabel('Data', 'FontSize', 30)
legend('Trajectory', 'Measurements', 'Exponential Mean', 'FontSize', 30)

%% 4) Second trajectory: Generate cyclic trajectory 𝑋𝑖 according to the equation 

A(1) = 1;
n_4 = 200;
T = 32;
omega = 2*pi/T;

sigma_w2= 0.08^2;
w = sqrt(sigma_w2).*randn(n_4,1);

for i = 2:n_4
    A(i) = A(i-1) + w(i);
end

x_sin = [];
for i = 1:n_4
    x_sin(i) = A(i) * sin(omega*i + 3);
end


%% 5) Generate measurements 𝑧𝑖 of the process 𝑋𝑖  

clear sigma_eta2
clear eta

sigma_eta2= 0.05;
eta = sqrt(sigma_eta2).*randn(n_4,1);

z_4 = [];
for i = 1:n_4
    z_4(i) = x_sin(i) + eta(i);
end

%% 6) Apply running mean with window size 𝑀=13 to measurements 𝑧𝑖. 

% Running mean
M_4 = 13;
x_hat_run_4 = movmean(z_4, M_4);

% plot of the results

figure(4)
plot(x_sin, 'r', 'LineWidth', 1.2)
hold on
plot(z_4, 'k', 'LineWidth', 1.2)
plot(x_hat_run_4, 'c', 'LineWidth', 1.2)
grid on; grid minor
xlabel('Steps', 'FontSize', 30)
ylabel('Data', 'FontSize', 30)
legend('Trajectory', 'Measurements', 'Running Mean', 'FontSize', 30)

%% Determine the period of oscillations for which running mean with given 
% for every group window size 𝑀

M_vect = [15:2:27];
x_hat_run = zeros(length(M_vect), 200);

for k = 1:length(M_vect)
      x_hat_run(k,:) = movmean(z_4, M_vect(k));
end

% plot
figure(5)
plot(x_sin, 'r', 'LineWidth', 1.2)
hold on
plot(z_4, 'k', 'LineWidth', 1.2)
plot(x_hat_run(end,:), 'c', 'LineWidth', 1.2)
grid on; grid minor
xlabel('Steps', 'FontSize', 30)
ylabel('Data', 'FontSize', 30)
legend('Trajectory', 'Measurements', 'Running Mean', 'FontSize', 30)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%
% TO DO:
%
% ask about deviation and variability indicators
% ask about M and alpha
% ask for the final point, not clear what to do!
%

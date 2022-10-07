%% Determining and removing drawbacks of exponential and running mean

close all; clear; clc;

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% I Part: Comparison of the traditional 13-month running mean with the forward-
%backward exponential smoothing for approximation of 11-year sunspot cycle

%% Download monthly mean sunspot numbe

data = load('data_group4.mat');
years = data.data(:,1); % year
m = data.data(:,2); % month
m_sun =  data.data(:,3); % monthly sunspot number
n = length(years);

%% 2) Make smoothing of monthly mean data by 13-month running mean

m_sun_mean = zeros(length(m_sun),1);
m_sun_mean(1:6) = 1/6*sum(m_sun(1:6));

m_sun_mean((length(m_sun_mean) - 5):length(m_sun_mean)) = 1/6 * sum(m_sun(length(m_sun)-5:length(m_sun)));

for i = 7:(n - 6)  
    m_sun_mean(i) = 1/24*(m_sun(i-6) + m_sun(i+6)) + 1/12*(m_sun(i-5) + m_sun(i-4) + m_sun(i-3) + m_sun(i-2) + m_sun(i-1) + ...
        m_sun(i) + m_sun(i+5) + m_sun(i+4) + m_sun(i+3) + m_sun(i+2) + m_sun(i+1));
end

 dev_13_step = sum((m_sun-m_sun_mean).^2);

 for i = 1:(n - 2)
    var_13_step(i) = (m_sun_mean(i + 2) - 2*m_sun_mean(i + 1) + m_sun_mean(i))^2;
 end
 var_13_step = sum(var_13_step);

%% 3) Make forward-backward exponential smoothing of monthly mean sunspot number. 
% Is there a smoothing constant 𝛼 that provides better results compared to 
% 13-month running mean according to deviation and variability indicators? 

alpha = [0.01, 0.02, 0.075, 0.1, 0.15, 0.2, 0.25];
x_forw = zeros(length(alpha), n);
x_forw(:,1) = m_sun(1);

x_back = zeros(length(alpha), n);
x_back(:,n) = x_forw(end);

for j = 1:length(alpha)
    for i = 2:n
        x_forw(j,i) = x_forw(j,i - 1) + alpha(j)*(m_sun(i) - x_forw(j,i - 1));
    end
    for i = (n - 1):-1:1
        x_back(j,i) = x_back(j,i + 1) + alpha(j)*(x_forw(j,i) - x_back(j,i + 1));
    end
end

% Variability indicators for different alpha
var_forw_sum = [];

for j = 1:length(alpha)
    for i = 1:(n - 2)
        var_forw(j,i) = (x_forw(j,i + 2) - 2*x_forw(j,i + 1) + x_forw(j,i))^2;
    end
    var_forw_sum(j) = sum(var_forw(j,:));
end



% Deviation indicators for different alpha
dev_forw = zeros(length(alpha), n);
dev_forw_sum = [];

for j = 1:length(alpha)
    for i = 1:n
        dev_forw(j,i) = (m_sun(i) - x_forw(j,i))^2;
    end
    dev_forw_sum(j) = sum(dev_forw(j,:));
end

figure(1)
plot(m_sun, 'c', 'LineWidth', 1.2)
hold on
plot(m_sun_mean, 'k', 'LineWidth', 1.2)
% plot(x_back(1,:), 'LineWidth', 1.2)
% plot(x_back(2,:), 'LineWidth', 1.2)
% plot(x_back(3,:), 'LineWidth', 1.2)
% plot(x_back(4,:), 'LineWidth', 1.2)
plot(x_back(5,:), 'LineWidth', 1.2)
grid on; grid minor
xlabel('Steps', 'FontSize', 30)
ylabel('Data', 'FontSize', 30)
% legend('m_sun', '13 smothing','a1', 'a2', 'a3', 'a4', 'a5', 'FontSize', 30)

%% II Part: 3d surface filtration using forward-backward smoothing 
%% 1) Download surface data

tr = load('true_surface.mat');
no = load('noisy_surface.mat');

true = tr.true_surface;
noise = no.noisy_surface;
x = [1:1:length(true)];

figure(1)
mesh(x,x,true)
colormap jet
xlabel('Y')
ylabel('X')
zlabel('Z')

figure(2)
mesh(x,x,noise)
colormap jet
xlabel('Y')
ylabel('X')
zlabel('Z')

%% 3) Determine the variance of deviation of noisy surface from the true one. 

var = sum((reshape((true-noise),...
    [1,length(true)*length(true)])).^2)/(length(true)*length(true)-1);

%% 4) 
alpha = 0.335;





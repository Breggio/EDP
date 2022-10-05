%% The relationship between solar radio flux F10.7 and sunspot number 

close all; clear all; clc;

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 1) Download monthly mean sunspot number and solar radio flux F10.7cm

data = load('data_group8.mat');
years = data.data_group8(:,1); % year
m = data.data_group8(:,2); % month
m_flux =  data.data_group8(:,3); % monthly solar radio flux at 10.7 cm 
m_sun =  data.data_group8(:,4); % monthly sunspot number
index = [1:length(years)];

%% 2) Plot the monthly mean sunspot number and solar radio flux F10.7 cm 
%     for visual representation. 

figure(1)
plot(index, m_flux, 'c', 'LineWidth', 1.2)
hold on
plot(index, m_sun, 'm', 'LineWidth', 1.2)
grid on; grid minor
xlabel('Years', 'FontSize', 30)
ylabel('Solar activity indicator', 'FontSize', 30)
legend('Solar radio flux F10.7 cm', 'Sunspot number', 'FontSize', 30) 
xlim([0 length(index)])
xticks([linspace(1,length(index), 12)])
xticklabels({'1953','1957','1961','1965','1969','1973','1977','1981',...
     '1985', '1989', '1993', '1997'})

% hold on
% figure(1)
% 
% yyaxis left
% plot(index, m_flux, 'k', 'LineWidth', 1.2)
% grid on; grid minor
% xlabel('Years', 'FontSize', 30)
% ylabel('Solar radio flux F10.7 cm', 'FontSize', 30)
% 
% yyaxis right
% plot(index, m_sun, 'm', 'LineWidth', 1.2)
% grid on; grid minor
% ylabel('Sunspot number', 'FontSize', 30)
% 
% legend('Solar radio flux F10.7 cm', 'Sunspot number', 'FontSize', 30) 
% xlim([1 526])
% hold off
%add the title in latex

%% 3) Make scatter plot between monthly mean sunspot number and solar radio flux F10.7 cm
figure(2)
scatter(m_sun, m_flux, 'b', 'LineWidth', 1.2) 
grid on; grid minor
xlabel('Sunspot number', 'FontSize', 30)
ylabel('Solar radio flux F10.7 cm', 'FontSize', 30)
xlim([0 length(index)])
xticks([linspace(1,length(index), 12)])
xticklabels({'1953','1957','1961','1965','1969','1973','1977','1981',...
     '1985', '1989', '1993', '1997'})
% we see a LINEAR correlation betweeen the number of sunspots and the solar radio flux 

%% 4) Make smoothing of monthly mean data (sunspot number and solar radio flux F10.7) 
%     by 13-month running mean.

m_sun_mean = zeros(length(m_sun),1);
m_flux_mean = zeros(length(m_flux),1);

m_sun_mean(1:6) = 1/6*sum(m_sun(1:6));
m_flux_mean(1:6) = 1/6*sum(m_flux(1:6));

m_sun_mean((length(m_sun_mean) - 5):length(m_sun_mean)) = 1/6 * sum(m_sun(length(m_sun)-5:length(m_sun)));
m_flux_mean((length(m_flux_mean) - 5):length(m_flux_mean)) = 1/6 * sum(m_flux(length(m_flux)-5:length(m_flux)));

for i = 7:(length(index) - 6)  
    m_sun_mean(i) = 1/24*(m_sun(i-6) + m_sun(i+6)) + 1/12*(m_sun(i-5) + m_sun(i-4) + m_sun(i-3) + m_sun(i-2) + m_sun(i-1) + ...
        m_sun(i) + m_sun(i+5) + m_sun(i+4) + m_sun(i+3) + m_sun(i+2) + m_sun(i+1));
    m_flux_mean(i) = 1/24*(m_flux(i-6) + m_flux(i+6)) + 1/12*(m_flux(i-5) + m_flux(i-4) + m_flux(i-3) + m_flux(i-2) + m_flux(i-1) + ...
        m_flux(i) + m_flux(i+5) + m_flux(i+4) + m_flux(i+3) + m_flux(i+2) + m_flux(i+1));
end

figure(3)
plot(index, m_flux_mean, 'k', 'LineWidth', 1.2)
hold on
plot(index, m_sun_mean, 'm', 'LineWidth', 1.2)
grid on; grid minor
xlabel('Years', 'FontSize', 30)
ylabel('Solar activity indicator', 'FontSize', 30)
legend('Solar radio flux F10.7 cm', 'Sunspot number', 'FontSize', 30) 
xlim([0 length(index)])
xticks([linspace(1,length(index), 12)])
xticklabels({'1953','1957','1961','1965','1969','1973','1977','1981',...
     '1985', '1989', '1993', '1997'})

% figure(3)
% hold on
% 
% yyaxis left
% plot(index, m_flux_mean, 'k', 'LineWidth', 1.2)
% grid on; grid minor
% xlabel('Years', 'FontSize', 30)
% ylabel('Solar radio flux F10.7 cm', 'FontSize', 30)
% 
% yyaxis right
% plot(index, m_sun_mean, 'm', 'LineWidth', 1.2)
% grid on; grid minor
% ylabel('Sunspot number', 'FontSize', 30)
% 
% legend('Solar radio flux F10.7 cm', 'Sunspot number', 'FontSize', 30) 
% xlim([1 526])
% hold off
%add the title in latex

%% 6-7) Construction of multi-dimensional linear regression

% m_flux_mean is the vector F of dependent variables (Vector of regressands)

% Matrix of independent variables (Vector of regressors)
R = ones(length(m_sun_mean),4); 
R(:,2) = m_sun_mean;
R(:,3) = m_sun_mean.^2;
R(:,4) = m_sun_mean.^3;

R_raw = ones(length(m_sun),4); 
R_raw(:,2) = m_sun;
R_raw(:,3) = m_sun.^2;
R_raw(:,4) = m_sun.^3;

%% 8) Determine vector of coefficients by LSM

% Vector of coefficients
beta_hat = (transpose(R)*R)^(-1)*transpose(R)*m_flux_mean;
beta_hat_raw = (transpose(R_raw)*R_raw)^(-1)*transpose(R_raw)*m_flux;

%% 9) Reconstruct solar radio flux at 10.7 cm on the basis of sunspot number

m_flux_mean_new = R * beta_hat;
m_flux_raw_new = R_raw * beta_hat_raw;

%% 10) Determine the variance of estimation error of solar radio flux

diff = [];
diff_raw = [];

for i = 1:length(m_flux_mean)

    diff = [diff; (m_flux_mean(i)-m_flux_mean_new(i))^2];
    diff_raw = [diff_raw; (m_flux(i)-m_flux_raw_new(i))^2];

end

var = 1/(length(m_flux_mean_new) - 1) * sum(diff)
var_raw = 1/(length(m_flux_raw_new) - 1) * sum(diff_raw)

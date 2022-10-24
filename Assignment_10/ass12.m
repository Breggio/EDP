%% Joint Assimilation of Navigation Datan Coming from Different Sources

% Written by Irina Yareshko and Luca Breggion, Skoltech 2022

close all 
clear
clc

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 

T=2; 
N=500;

InititalState = [1000;100;1000;100]; %[x0; Vx0; y0; Vy0]

sigma_a = 0.3;          % variance of acceleration noise
sigma_D = 50 ;           % variance of range noise of measurements
sigma_b = 0.004;        % variance of azimuth noise of measurements
sigma_b_add = 0.001;    %

% Creating of arrays for running Kalman filter M times

M=500;                            %number of runs
err_range_filt  = zeros(M,N); %Filtraion error of range
err_range_pred  = zeros(M,N); %Predicted error of range
err_azimut_filt = zeros(M,N); %Filtraion error of azimut
err_azimut_pred = zeros(M,N); %Predicted error of azimut

for i=1:M
    %Generation of deterministic trajectory and its measurements
    [polar,cart,z_p] = ...
        true_traj(N,T,sigma_D,sigma_b,InititalState, sigma_a, sigma_b_add);

    %Applying or Kalman filter for measurements
    [Z_filtered,Z_forecast, range_fe, azimuth_fe] = Kalman_extended(z_p,...
        T,sigma_D,sigma_b, sigma_a, sigma_b_add);
    
    err_range_filt(i,:) = (polar(1,:) - range_fe(1,:)).^2;
    err_range_pred(i,:) = (polar(1,:) - range_fe(2,:)).^2;
    err_azimut_filt(i,:)= (polar(2,:) - azimuth_fe(1,:)).^2;
    err_azimut_pred(i,:)= (polar(2,:) - azimuth_fe(2,:)).^2;
end

err_range_filt =sqrt(1/(M-1)*sum(err_range_filt));
err_range_pred =sqrt(1/(M-1)*sum(err_range_pred));
err_azimut_filt = sqrt(1/(M-1)*sum(err_azimut_filt));
err_azimut_pred = sqrt(1/(M-1)*sum(err_azimut_pred));

%% Point 6

% True trajectory, measurements and filtered&extrapolated
figure(1)
polarplot(polar(2,:),polar(1,:), 'm',...
    z_p(2,:),z_p(1,:),'c.',...
    azimuth_fe(1,:), range_fe(1,:), 'k')
legend('True trajectory', 'Measurements', 'Filtered', 'FontSize', 30)
grid on; grid minor

%% Point 7

figure(2)
plot(4:N,err_range_filt(4:N),'m',...
     4:N,err_range_pred(4:N),'c',...
     1:N,sigma_D*ones(1,N), 'black');
legend('Filtration error','Extrapolation error','$\sigma_D$','interpreter', 'latex');
%title('Errors of range');
xlabel('Step')
ylabel('Errors')
grid on; grid minor

%% Point 7

figure(3)
plot(4:N,err_azimut_filt(4:N),'m', ...
     4:N,err_azimut_pred(4:N),'c',...
     1:N,sigma_b*ones(1,N), 'black');
legend('Filtration error','Extrapolation error','$\sigma_\beta$', 'interpreter', 'latex');
%title('Errors of azimiuth');
xlabel('Step')
ylabel('Errors')
grid on; grid minor

%% Point 7

figure(4)
plot(cart(1,:),cart(3,:), 'm',...
    Z_filtered(1,:),Z_filtered(3,:),'c')
%title(True and filtered trajectory')
legend('True', 'Filtered')
grid on; grid minor
xlabel('Step')
ylabel('Errors')

%% Point 8

figure(7)
plot(1:N, polar(1,:),'m',  1:N, z_p(1,:), 'c',1:N ,range_fe(1,:), 'k')
%title('Range')
legend('True', 'Measurements', 'Filtered')
grid on; grid minor
xlabel('Step')
ylabel('Data')

%% Point 8

figure(8)
plot(1:N, polar(2,1:N),'m', 4:N, z_p(2,4:N), 'c', 1:N, azimuth_fe(1,:), 'k')
%title('Azimuth')
legend('True', 'Measurements', 'Filtered')
grid on; grid minor
xlabel('Step')
ylabel('Data')

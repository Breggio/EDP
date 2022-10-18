% Written by Irina Yareshko and Luca Breggion, Skoltech 2022

close all 
clear
clc

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%Initial data
T = 2; N = 26;
InititalState = [13500/sqrt(2); -50; 13500/sqrt(2); -45]; 
sigma_D = 20;      %variance of range noise of measurements
sigma_beta = 0.02; %variance of azimuth noise of measurements

%Creating of arrays for running Kalman filter M times
M = 500;                          %number of runs
Err_range_filtered  = zeros(M,N); %Filtraion error of range
Err_range_forecast  = zeros(M,N); %Predicted error of range
Err_azimut_filtered = zeros(M,N); %Filtraion error of azimut
Err_azimut_forecast = zeros(M,N); %Predicted error of azimut
Condition_nums      = zeros(M,N); %array of average condition number
K_matr_values       = zeros(M,N); 

Mean_Z_x     = zeros(1,N);
Mean_azimuth = zeros(1,N);

Mean_range = zeros(1,N);
for i=1:M
    %Generation of deterministic trajectory and its measurements
    [TruePolar, TrueCart, Z_c, Z_p] = ...
        Generation_true_determ(N, T, sigma_D, sigma_beta, InititalState);
    %Applying or Kalman filter for measurements
    [Z_filtered, P, P_pred, ...
        range_fe, azimuth_fe, CondMatr, K_matr] = ...
        Kalman(Z_c, Z_p, T, sigma_D, sigma_beta);
    
    Err_range_filtered(i,:)  = (TruePolar(1,:) - range_fe(1,:)).^2;
    Err_range_forecast(i,:)  = (TruePolar(1,:) - range_fe(2,:)).^2;
    Err_azimut_filtered(i,:) = (TruePolar(2,:) - azimuth_fe(1,:)).^2;
    Err_azimut_forecast(i,:) = (TruePolar(2,:) - azimuth_fe(2,:)).^2;
    
    Condition_nums(i,:) = CondMatr;
    K_matr_values(i,:) = K_matr(1,1,:);
    
    Mean_Z_x = Mean_Z_x + Z_filtered(1,2:N+1);
    Mean_azimuth = Mean_azimuth + azimuth_fe(1,:);
end
Mean_Z_x = Mean_Z_x/M;
Mean_azimuth = Mean_azimuth/M;

FinalErr_range_filtered  = sqrt(1/(M - 1)*sum(Err_range_filtered));
FinalErr_range_forecast  = sqrt(1/(M - 1)*sum(Err_range_forecast));
FinalErr_azimut_filtered = sqrt(1/(M - 1)*sum(Err_azimut_filtered));
FinalErr_azimut_forecast = sqrt(1/(M - 1)*sum(Err_azimut_forecast));
Final_CN = sum(Condition_nums)/M;
Final_K_matr = sum(K_matr_values)/M;

%% Point 2 (Generated motion in polar coordinate system)
figure(1)
polarplot(TruePolar(2,:), TruePolar(1,:), 'LineWidth', 1.2)
grid on; grid minor
legend('True motion', 'FontSize', 20);
title('Object moves uniformly', 'FontSize', 20);

%% Point 10 (Errors of extrapolation and filtration estimates of range and azimuth)
figure(2)
plot(3:N,FinalErr_range_filtered(3:N),...
     3:N,FinalErr_range_forecast(3:N),...
     3:N,sigma_D*ones(1, N-2), 'black', 'LineWidth', 1.2);
grid on; grid minor
legend('True filtration error','True extrapolation error','\sigma_D', 'FontSize', 20);
title('a) Errors of range', 'FontSize', 20);
xlabel('Step', 'FontSize', 20)
ylabel('Errors', 'FontSize', 20)

figure(3)
plot(3:N,FinalErr_azimut_filtered(3:N), ...
     3:N,FinalErr_azimut_forecast(3:N), ...
     3:N,sigma_beta*ones(1, N-2), 'black', 'LineWidth', 1.2);
grid on; grid minor
legend('True filtration error','True extrapolation error','\sigma_\beta', 'FontSize', 20);
title('b) Errors of azimuth', 'FontSize', 20);
xlabel('Step', 'FontSize', 20)
ylabel('Errors', 'FontSize', 20)

%% Point 12 (Building of plot of dependence of coordinate x on azimuth b)
figure(4)
plot(Mean_azimuth(1,:), Mean_Z_x(1,:),'r',TruePolar(2,:), TrueCart(1,:), 'b--', 'LineWidth', 1.2);
grid on; grid minor
legend('x = f(\beta) filtered values', 'x = g(\beta) true data', 'FontSize', 20);
title('Dependance of coordinate x on azimuth \beta', 'FontSize', 20);
xlabel('Azimuth', 'FontSize', 20)
ylabel('Coordinate x', 'FontSize', 20)

%% Point 12 (Building of plot of dinamics of condotion number)
figure(5)
plot(1:N, Final_CN, 'LineWidth', 1.2);
grid on; grid minor
legend('Condition number', 'FontSize', 20);
title('Dinamics of condition number', 'FontSize', 20);
xlabel('Step', 'FontSize', 20)
ylabel('Condition number', 'FontSize', 20)
%% Point 13 (Building of plot of dinamics of filter gain)
figure(6)
plot(1:N, Final_K_matr, 'LineWidth', 1.2);
grid on; grid minor
legend('Filter gain','FontSize', 20);
title('Dinamics of filter gain', 'FontSize', 20);
xlabel('Step','FontSize', 20)
ylabel('K','FontSize', 20)

%% Point 14 (Quite close distance from an observer)
InititalState = [3500/sqrt(2); -50; 3500/sqrt(2); -45];
for i=1:M
    %Generation of deterministic trajectory and its measurements
    [TruePolar, TrueCart, Z_c, Z_p] = ...
        Generation_true_determ(N, T, sigma_D, sigma_beta, InititalState);
    %Applying or Kalman filter for measurements
    [Z_filtered, P, P_pred, ...
        range_fe, azimuth_fe, CondMatr, K_matr] = ...
        Kalman(Z_c, Z_p, T, sigma_D, sigma_beta);
    
    Err_range_filtered(i,:)  = (TruePolar(1,:) - range_fe(1,:)).^2;
    Err_range_forecast(i,:)  = (TruePolar(1,:) - range_fe(2,:)).^2;
    Err_azimut_filtered(i,:) = (TruePolar(2,:) - azimuth_fe(1,:)).^2;
    Err_azimut_forecast(i,:) = (TruePolar(2,:) - azimuth_fe(2,:)).^2;
    
    Condition_nums(i,:) = CondMatr;
    K_matr_values(i,:) = K_matr(1,1,:);
    
    Mean_Z_x = Mean_Z_x + Z_filtered(1,2:N+1);
    Mean_azimuth = Mean_azimuth + azimuth_fe(1,:);
    Mean_range = Mean_range + range_fe(1,:);
end
Mean_Z_x = Mean_Z_x/M;
Mean_azimuth = Mean_azimuth/M;
Mean_range = Mean_range/M;

FinalErr_range_filtered  = sqrt(1/(M-1)*sum(Err_range_filtered));
FinalErr_range_forecast  = sqrt(1/(M-1)*sum(Err_range_forecast));
FinalErr_azimut_filtered = sqrt(1/(M-1)*sum(Err_azimut_filtered));
FinalErr_azimut_forecast = sqrt(1/(M-1)*sum(Err_azimut_forecast));
Final_CN = sum(Condition_nums)/M;

% Generated motion in polar coordinate system (Quite close)
figure(7)
polarplot(TruePolar(2,:), TruePolar(1,:), 'LineWidth', 1.2)
grid on; grid minor
legend('True motion');
title('Object moves uniformly (Quite close)');

%% Point 15 (Errors of extrapolation and filtration estimates of (Quite close))
figure(8)
subplot(1,2,1);
plot(3:N, FinalErr_range_filtered(3:N),...
     3:N, FinalErr_range_forecast(3:N),...
     3:N, sigma_D*ones(1, N-2), 'black', 'LineWidth', 1.2);
grid on; grid minor
legend('True filtration error','True extrapolation error','\sigma_D');
title('a) Errors of range (Quite close)');
xlabel('Step')
ylabel('Errors')

figure(9)
plot(3:N, FinalErr_azimut_filtered(3:N),...
     3:N, FinalErr_azimut_forecast(3:N),...
     3:N, sigma_beta*ones(1, N-2), 'black');
legend('True filtration error','True extrapolation error','\sigma_\beta');
title('b) Errors of azimuth (Quite close)');
xlabel('Step')
ylabel('Errors')

%% Point 16 (Dependence of coordinate x on azimuth b (Quite close))
figure(10)
plot(Mean_azimuth(1,:), Mean_Z_x(1,:),'r',TruePolar(2,:), TrueCart(1,:), 'b--', 'LineWidth', 1.2);
grid on; grid minor
legend('x = f(\beta) filtered values', 'x = g(\beta) true data','FontSize', 20);
title('Dependance of coordinate x on azimuth \beta (Quite close)','FontSize', 20);
xlabel('Azimuth','FontSize', 20)
ylabel('Coordinate x','FontSize', 20)

%% Point 17 (Building of plot of dinamics of condotion number (Quite close))
figure(11)
plot(1:N,Final_CN, 'LineWidth', 1.2);
grid on; grid minor
legend('Condition number','FontSize', 20);
title('Dinamics of condition number (Quite close)','FontSize', 20);
xlabel('Step','FontSize', 20)
ylabel('Condition number','FontSize', 20)
%% Point 19 quite close distance from an observer and other values of variances
sigma_D=50; %variance of range noise of measurements
sigma_beta=0.0015; %variance of azimuth noise of measurements
for i=1:M
    %Generation of deterministic trajectory and its measurements
    [TruePolar,TrueCart,Z_c,Z_p] = ...
        Generation_true_determ(N,T,sigma_D,sigma_beta,InititalState);
    %Applying or Kalman filter for measurements
    [Z_filtered, P, P_pred, ...
        range_fe, azimuth_fe, CondMatr, K_matr] = ...
        Kalman(Z_c,Z_p,T,sigma_D,sigma_beta);
    
    Err_range_filtered(i,:) = (TruePolar(1,:) - range_fe(1,:)).^2;
    Err_range_forecast(i,:) = (TruePolar(1,:) - range_fe(2,:)).^2;
    Err_azimut_filtered(i,:)= (TruePolar(2,:) - azimuth_fe(1,:)).^2;
    Err_azimut_forecast(i,:)= (TruePolar(2,:) - azimuth_fe(2,:)).^2;
    
    Condition_nums(i,:) = CondMatr;
    K_matr_values(i,:) = K_matr(1,1,:);
    
    Mean_Z_x = Mean_Z_x + Z_filtered(1,2:N+1);
    Mean_azimuth = Mean_azimuth + azimuth_fe(1,:);
end
Mean_Z_x = Mean_Z_x/M;
Mean_azimuth = Mean_azimuth/M;

FinalErr_range_filtered =sqrt(1/(M-1)*sum(Err_range_filtered));
FinalErr_range_forecast =sqrt(1/(M-1)*sum(Err_range_forecast));
FinalErr_azimut_filtered = sqrt(1/(M-1)*sum(Err_azimut_filtered));
FinalErr_azimut_forecast = sqrt(1/(M-1)*sum(Err_azimut_forecast));
Final_CN=sum(Condition_nums)/M;
Final_K_matr=sum(K_matr_values)/M;

% Generated motion in polar coordinate system (Quite close, other variances)
figure(12)
polarplot(TruePolar(2,:),TruePolar(1,:), 'LineWidth', 1.2)
grid on; grid minor
legend('True motion','FontSize', 20);
title('Object moves uniformly (Quite close, other variances)','FontSize', 20);

%Errors of extrapolation and filtration estimates of (Quite close, other variances)
figure(13)
plot(3:N,FinalErr_range_filtered(3:N),...
     3:N,FinalErr_range_forecast(3:N),...
     3:N,sigma_D*ones(1,N-2), 'black', 'LineWidth', 1.2);
grid on; grid minor
legend('True filtration error','True extrapolation error','\sigma_D','FontSize', 20);
title('a) Errors of range (Quite close, other variances)','FontSize', 20);
xlabel('Step','FontSize', 20)
ylabel('Errors','FontSize', 20)

figure(14)
plot(3:N,FinalErr_azimut_filtered(3:N),...
     3:N,FinalErr_azimut_forecast(3:N),...
     3:N,sigma_beta*ones(1,N-2), 'black', 'LineWidth', 1.2);
grid on; grid minor
legend('True filtration error','True extrapolation error','\sigma_\beta','FontSize', 20);
title('b) Errors of azimuth (Quite close, other variances)','FontSize', 20);
xlabel('Step','FontSize', 20)
ylabel('Errors','FontSize', 20)

%Dependence of coordinate x on azimuth b (Quite close, other variances)
figure(15)
plot(Mean_azimuth(1,:), Mean_Z_x(1,:),'r',TruePolar(2,:), TrueCart(1,:), 'b--', 'LineWidth', 1.2);
grid on; grid minor
legend('x = f(\beta) filtered values', 'x = g(\beta) true data','FontSize', 20);
title('Dependance of coordinate x on azimuth \beta (Quite close, other variances)','FontSize', 20);
xlabel('Azimuth','FontSize', 20)
ylabel('Coordinate x','FontSize', 20)

% Building of plot of dinamics of condotion number (Quite close, other variances)
figure(16)
plot(1:N, Final_CN, 'LineWidth', 1.2);
grid on; grid minor
legend('Condition number','FontSize', 20);
title('Dinamics of condition number (Quite close, other variances)','FontSize', 20);
xlabel('Step','FontSize', 20)
ylabel('Condition number','FontSize', 20)
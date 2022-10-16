%% Development of optimal smoothing to increase the estimation accuracy 

% Written by Irina Yareshko and Luca Breggion, Skoltech 2022

close all 
clear
clc

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Initial data
N = 200; %observation interval
sigma2_n = 20^2; %variance of noise
sigma2_a = 0.2^2; %variance of acceleraration
T=1; %period of step

M=500; %number of runs
X0 = [2; 0]; P0 = [10000 0; 0 10000];
SqrError_X = zeros(M, 1, N - 2); %Squared error of coordinate
SqrError_V = zeros(M, 1, N - 2); %Squared error of velocity
for i=1:M
    %Generation of data
    [X, z, V] = data_gen(N,T,sigma2_n,sigma2_a);
    %Applying of Kalman filter
    [Z_f, X_f, P_pred, P] = kalman(z,  sigma2_a, X0, P0);
    [Z_smoothed, P_smooth] = kalman_sm(Z_f, P_pred, P, T); %smoothing of filtered data
    SqrError_X(i,:) = ( X(3:N) - Z_smoothed(1,4:N+1) ).^2;
    SqrError_V(i,:) = ( V(3:N) - Z_smoothed(2,4:N+1) ).^2;
end

Final_ErrSmoothed_X = sqrt( 1/(M-1) * sum(SqrError_X) ); %true estimation error of cordinate
Final_ErrSmoothed_V = sqrt( 1/(M-1) * sum(SqrError_V) ); %true estimation error of velocity

SmoothError_X = sqrt(P_smooth(1,1,4:N + 1)); %smoothed error of cordinate
SmoothError_V = sqrt(P_smooth(2,2,4:N + 1)); %smoothed error of velocity

t = 1:N;
figure(1)
plot(t, X, 'g', t, z, 'm-', t, Z_f(1,2:N + 1), 'k', t, Z_smoothed(1,2:N + 1), 'c', 'LineWidth', 1.2)
grid on; grid minor
legend('True Data', 'Measurements', 'Filtered Data', 'Smoothed Data', 'FontSize', 20);
xlabel('Step', 'FontSize', 30);
ylabel('Data', 'FontSize', 30);

figure(2)
plot(3:N, SmoothError_X(1,:), 'm', 3:N, Final_ErrSmoothed_X(1,:), 'k', 'LineWidth', 1.2)
grid on; grid minor
legend('Smoothing Algorithm Error', 'True Estimation Error', 'FontSize', 20);
xlabel('Step', 'FontSize', 30);
ylabel('Error', 'FontSize', 30)

figure(3)
plot(3:N, SmoothError_V(1,:), 'm', 3:N, Final_ErrSmoothed_V(1,:), 'k', 'LineWidth', 1.2)
grid on; grid minor
legend('Smoothing Algorithm Error', 'True Estimation Error', 'FontSize', 20)
xlabel('Step', 'FontSize', 30);
ylabel('Error', 'FontSize', 30)

FilteredError_X = sqrt(P(1,1,4:N + 1)); %filtered error of cordinate
FilteredError_V = sqrt(P(2,2,4:N + 1)); %filtered error of velocity

figure(4)
plot(3:N, SmoothError_X(1,:), 'm', 3:N, FilteredError_X(1,:), 'k', 'LineWidth', 1.2)
grid on; grid minor
legend('Smoothing Algorithm Error', 'Filteration Error', 'FontSize', 20);
xlabel('Step', 'FontSize', 30);
ylabel('Error', 'FontSize', 30)

figure(5)
plot(3:N, SmoothError_V(1,:), 'm', 3:N, FilteredError_V(1,:), 'k', 'LineWidth', 1.2)
grid on; grid minor
legend('Smoothing Algorithm Error', 'Filteration Error', 'FontSize', 20);
xlabel('Step', 'FontSize', 30);
ylabel('Error', 'FontSize', 30)


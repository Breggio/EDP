%% Tracking and forecasting in conditions of measurement gaps 

% Written by Irina Yareshko and Luca Breggion, Skoltech 2022

close all 
clear
clc

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 1-2) Develop Kalman filter to track moving object under conditions of gap 𝑃=0.2

sigma2_a = 0.2^2;
sigma2_n = 20^2;

%Initial data
N = 200; % observation interval
t = 1:N; %$ observation interval vector
T = 1;
Prob = 0.2;

[x, z, x_Kalman] = tracking(N, T, sigma2_a, sigma2_n, Prob);

figure(1)
plot(t, x, 'g', t, z, 'm-', t, x_Kalman,'k', 'LineWidth',1.5)
grid on; grid minor
legend('True Data', 'Measurements with P = 0.2', 'Kalman Filter with P = 0.2', 'FontSize', 20, 'location', 'best')
xlabel('Step', 'FontSize', 30)
ylabel('Data', 'FontSize', 30)

%% 3) Determine filtered and extrapolated errors of estimation
% (1 step and 7 steps ahead) over 500 runs of filter. 
% Compare them with true estimation errors. 

%Initial data
m = 7;            % number of steps ahead
M = 500;          % number of runs
N = 200;          % observation interval
sigma2_n = 20^2;  % variance of noise
sigma2_a = 0.2^2; % variance of acceleraration
T = 1;            % period of step

[Final_Error_X, Final_Error_X_f, Final_Error_X_f7] = errors(m, M, N,sigma2_n, sigma2_a, T, x);

t1 = 7:N;

figure(2)
plot(t, Final_Error_X(1:N),'m',t, Final_Error_X_f, 'c', t1, Final_Error_X_f7,'k', 'LineWidth', 1.2);
grid on; grid minor
xlabel('Step', 'FontSize', 30)
ylabel('Errors', 'FontSize', 30)
legend('True estimation error','True error of 1-step prediction with P=0.2', ...
    'True error of 7-step prediction with P=0.2', 'FontSize', 30)

%% 4) Analyze the decrease of estimation accuracy in conditions of measurement gaps. 

P1 = 0.3;
P2 = 0.5;
P3 = 0.7;

% Tracking with different P
[x1, z1, x_Kalman1] = tracking(N, T, sigma2_a, sigma2_n, P1);
[x2, z2, x_Kalman2] = tracking(N, T, sigma2_a, sigma2_n, P2);
[x3, z3, x_Kalman3] = tracking(N, T, sigma2_a, sigma2_n, P3);

% Errors
[Final_Error_X1, Final_Error_X_f1, Final_Error_X_f71] = errors(m, M, N,sigma2_n, sigma2_a, T, x1);
[Final_Error_X2, Final_Error_X_f2, Final_Error_X_f72] = errors(m, M, N,sigma2_n, sigma2_a, T, x2);
[Final_Error_X3, Final_Error_X_f3, Final_Error_X_f73] = errors(m, M, N,sigma2_n, sigma2_a, T, x3);

%% Plots

figure(3)
plot(t, x1, 'g', t, z1, 'm-', t, x_Kalman1,'k', 'LineWidth',1.5)
grid on; grid minor
legend('True Data', 'Measurements with P = 0.3', 'Kalman Filter with P = 0.3', 'FontSize', 20, 'location', 'best')
xlabel('Step', 'FontSize', 30)
ylabel('Data', 'FontSize', 30)

figure(4)
plot(t, x2, 'g', t, z2, 'm-', t, x_Kalman2,'k', 'LineWidth',1.5)
grid on; grid minor
legend('True Data', 'Measurements with P = 0.5', 'Kalman Filter with P = 0.5', 'FontSize', 20, 'location', 'best')
xlabel('Step', 'FontSize', 30)
ylabel('Data', 'FontSize', 30)

figure(5)
plot(t, x3, 'g', t, z3, 'm-', t, x_Kalman3,'k', 'LineWidth',1.5)
grid on; grid minor
legend('True Data', 'Measurements with P = 0.7', 'Kalman Filter with P = 0.7', 'FontSize', 20, 'location', 'best')
xlabel('Step', 'FontSize', 30)
ylabel('Data', 'FontSize', 30)

figure(6)
plot(t, Final_Error_X1(1:N),'m',t, Final_Error_X_f1, 'c', t1, Final_Error_X_f71,'k', 'LineWidth', 1.2);
grid on; grid minor
xlabel('Step', 'FontSize', 30)
ylabel('Errors', 'FontSize', 30)
legend('True estimation error','True error of 1-step prediction with P=0.3', ...
    'True error of 7-step prediction with P=0.3', 'FontSize', 30)

figure(7)
plot(t, Final_Error_X2(1:N),'m',t, Final_Error_X_f2, 'c', t1, Final_Error_X_f72,'k', 'LineWidth', 1.2);
grid on; grid minor
xlabel('Step', 'FontSize', 30)
ylabel('Errors', 'FontSize', 30)
legend('True estimation error','True error of 1-step prediction with P=0.5', ...
    'True error of 7-step prediction with P=0.5', 'FontSize', 30)

figure(8)
plot(t, Final_Error_X3(1:N),'m',t, Final_Error_X_f3, 'c', t1, Final_Error_X_f73,'k', 'LineWidth', 1.2);
grid on; grid minor
xlabel('Step', 'FontSize', 30)
ylabel('Errors', 'FontSize', 30)
legend('True estimation error','True error of 1-step prediction with P=0.7', ...
    'True error of 7-step prediction with P=0.7', 'FontSize', 30)

figure(9)
plot(t, Final_Error_X1(1:N),'m', t, Final_Error_X2(1:N),'c', t, Final_Error_X3(1:N),'k', 'LineWidth',1.5)
grid on; grid minor
xlabel('Step', 'FontSize', 30)
ylabel('Errors', 'FontSize', 30)
legend('True estimation error with P = 0.3', 'True estimation error with P = 0.5',...
    'True estimation error with P = 0.7','FontSize', 30)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                FUNCTION                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, z, x_Kalman] = tracking(N, T, sigma2_a, sigma2_n, Prob)

eta = randn*sqrt(sigma2_n); % Random noise of measurments

x(1) = 5;
V(1) = 1;
Z(1) = x(1) + eta; %the first measurment

%State matrixes
phi = [1 T; 0 1];
G  = [T/2; T];    % input matrix
H  = [1 0];       % observation matrix

% Kalman filter 
X = zeros(2, 2, N + 1); %true data
P = zeros(2, 4, N + 1); %filtration error covariance matrix
K = zeros(2, N);
X_f = zeros(2, N);

% Initial conditions
X(:,:,1) = [2, 0; 0, 0];
Q = G*G.'*sigma2_a; 
R = sigma2_n;
m = 7; % steps ahead

for i=2:N

    a = randn*sqrt(sigma2_a); %normally distributed random acceleration
    eta = randn*sqrt(sigma2_n); %random noise of measurments
    V(i) = V(i - 1) + a*T;
    x(i) = x(i - 1) + V(i - 1)*T + a*T^2/2;

    csi = rand;

    if csi <= Prob
        z(i) = NaN;
        X(:,:,i) = X(:,:,i-1);
        P(:,:,i) = P(:,:,i-1);
        X_f = phi^m*X(:,1,i);
    
    elseif csi > Prob
        z(i) = x(i) + eta;
        X_f(:,i-1) = phi^7*X(:,1,i-1);
        X(:,2,i-1) = phi*X(:,1,i-1);
        P(:,3:4,i-1) = phi*P(:,1:2,i-1)*phi.'+Q;
        K(:,i-1) = P(:,3:4,i-1)*H.'*(H*P(:,3:4,i-1)*H.'+R) ^(-1);

        if isnan(z(i-1)) == 1 
            counter=i ;
            znew=z(i-1);
            
            while isnan(znew) == 1 
                znew=z ( counter-1); 
                counter=counter-1;
            end
        elseif isnan(z(i-1))== 0 
            znew=z(i-1);
        end
        
        X(:,1,i) = X(:,2,i-1)+K(:,i-1)*(znew-H*X(:,2, i-1) ) ;
        P(:,1:2,i) = (eye(2)-K(:,i-1)*H)*P(:,3:4,i-1);

    end
end

x_Kalman(1:N) = X(1,1,1:N);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                FUNCTION                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Final_Error_X, Final_Error_X_f, Final_Error_X_f7] = errors(m, M, N,sigma2_n, sigma2_a, T, x)

% Initialization
Error_X    = zeros(M, N);      % array of errors of filtered estimate
Error_X_f  = zeros(M, N);      % array of errors of forecasts
Error_X_f7 = zeros(M, N - 6);  % array of errors of forecasts
X_Kalman   = zeros(6, N, M);   % array of filtered data
%x          = zeros(2, N);

for i = 1:M
    [x, Z] = data_gen(N,T,sigma2_n,sigma2_a);
    X_Kalman(:,:,i) = Kalman_filter(Z,T,m,sigma2_n,sigma2_a); % Kalman filter
    Error_X(i,:)    = (x(1,:) - X_Kalman(1,:,i)).^2; % errors of filtered estimate
    Error_X_f(i,:)  = (x(1,:) - X_Kalman(2,:,i)).^2; % errors of forecasts 1-step
    Error_X_f7(i,:) = (x(1,7:N) - X_Kalman(3,1:N-6,i)).^2; % errors of forecasts 7-step
end

%Final average value of Error over M runs
Final_Error_X    = sqrt(1/(M-1)*sum(Error_X));
Final_Error_X_f  = sqrt(1/(M-1)*sum(Error_X_f));
Final_Error_X_f7 = sqrt(1/(M-1)*sum(Error_X_f7));

end




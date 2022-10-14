%% Tracking of a moving object which trajectory is disturbed by random acceleration 

close all; clear; clc;

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 1. Generate a true trajectory 𝑋𝑖 of an object motion disturbed by normally...
% distributed random acceleration 

N = 200;

x(1) = 5;
v(1) = 1;
t = 1;

sigma_a2 = 0.2^2; % variance of noise
sigma_eta2 = 20^2;

a = sqrt(sigma_a2).*randn(N,1);
eta = sqrt(sigma_eta2).*randn(N,1);

for i = 2:N
    v(i) = v(i - 1) + a(i - 1)*t;
    x(i) = x(i - 1) + v(i - 1)*t + (a(i - 1)*t^2) * 0.5;   
end

Z = [];

for i = 1:N
    Z(i) = x(i) + eta(i);
end

%% 2) To construct Kalman filter we need to present the system at state space 
% taking into account that only measurements of coordinate 𝑥𝑖 are available 

% Transition matrix
phi = [1, t; 0, 1];

% Input matrix
G = [t^2/2;t];

% Observation matrix
H = [1,0];

% Initial filtered estimate
X(:,:,1) = [2,0;0,0];

% Initial filtration error covariance matrix
P(:, 1:2, 1) = [1000, 0; 0, 1000];

% Covariance matrix of state noise
Q = G*transpose(G)*sigma_a2;

% Covariance matrix of measurements noise
R = sigma_eta2;

for i = 2:N+1
    X(:,2,i-1) = phi*X(:,1,i-1);
    X_f(:,i-1) = phi^7*X(:,1,i-1);
    P(:,3:4,i-1) = phi*P(:,1:2,i-1)*phi.'+Q;
    K(:,i-1)=P(:,3:4,i-1)*H.'*(H*P(:,3:4,i-1)*H.'+R)^(-1);
    X(:,1,i)=X(:,2,i-1)+K(:,i-1)*(Z(i-1)-H*X(:,2,i-1));
    P(:,1:2,i)=(eye(2)-K(:,i-1)*H)*P(:,3:4,i-1);
end

x_Kalman(1:N)=X(1,1,2:N+1);

t=1:N; %array of steps
figure
plot(t,x,'g',t,Z,'r',t,x_Kalman,'k')
legend('true data','measurments','filtered estimates of state vector');
grid on
xlabel('Step')
ylabel('Data')
title('Applying of backward exponential mean')

%% Tracking of a moving object which trajectory is disturbed by random acceleration 

close all; clear; clc;

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 1. Generate a true trajectory 𝑋𝑖 of an object motion disturbed by normally...
% distributed random acceleration 

n = 200;

x(1) = 5;
v(1) = 1;
t = 1;

sigma_a2 = 0.2^2; % variance of noise
sigma_eta2 = 20^2;

a = sqrt(sigma_a2).*randn(n,1);
eta = sqrt(sigma_eta2).*randn(n,1);

for i = 2:n
    v(i) = v(i - 1) + a(i - 1)*t;
    x(i) = x(i - 1) + v(i - 1)*t + (a(i - 1)*t^2) * 0.5;   
end

z = [];

for i = 1:n
    z(i) = x(i) + eta(i);
end

%% 2) To construct Kalman filter we need to present the system at state space 
% taking into account that only measurements of coordinate 𝑥𝑖 are available 

% Transition matrix
phi = [1, t; 0, 1];

% Input matrix
G = [t^2/2;t];

% Observation matrix
H = [1,0];

X_k = [x;v];

for i=1:n
    z_k(i) = H*X_k(:,i) + eta(i);
end
  
% Initial filtered estimate
X_k_2 = [2;0];

for i=2:n
    X_k_2(:,i) = phi*X_k(:,i-1) + G*a(i-1);
end

% Initial filtration error covariance matrix
P = [10000,0 ; 0, 10000];

% Covariance matrix of state noise
Q = G*transpose(G)*sigma_a2;

% Covariance matrix of measurements noise
R = sigma_eta2;




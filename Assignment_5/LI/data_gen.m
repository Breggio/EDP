function [X, Z] = data_gen(N,T,sigma2_n,sigma2_a)

X = zeros(1,N); % true data
V = zeros(1,N); % velocity
Z = zeros(1,N); % measurments

n = randn*sqrt(sigma2_n); %random noise of measurments

% Initial data
X(1) = 5;
V(1) = 1;
Z(1) = X(1) + n; % first measurment

for i = 2:N
    a = randn*sqrt(sigma2_a); % normally distributed random acceleration
    n = randn*sqrt(sigma2_n); % random noise of measurments
    V(i) = V(i-1) + a*T;
    X(i) = X(i-1) + V(i-1)*T + a*T^2/2;
    Z(i) = X(i) + n;
end
end
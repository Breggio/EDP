function [x, x_hat, z, alpha, sigma_eta2, sigma_w2] = fun(n, var_w, var_eta, incond)

a_1 = sqrt(var_w);
a_2 = sqrt(var_eta);

w = a_1.*randn(n,1);
eta = a_2.*randn(n,1);

x(1) = incond;

for i = 2:n
    x(i) = x(i-1) + w(i); % generated trajectory RWM
end

for i=1:n
    z(i) = x(i) + eta(i); % Generate measurements 𝑧𝑖 of the process 𝑋𝑖 
end

% 2) Identify  𝜎𝑤2  and  𝜎𝜂2  using  identification  method

E_v_sq_sum = [];
for i = 2:n
    E_v_sq_sum(i-1) = ( w(i) + eta(i) - eta(i-1) )^2;
end

E_v_sq = 1/(n-1) *sum(E_v_sq_sum);

E_rho_sq_sum = [];
for i=3:n
     E_rho_sq_sum(i-2) = ( w(i) + w(i-1) + eta(i) - eta(i-2) )^2;
end

E_rho_sq = 1/(n-2) *sum(E_rho_sq_sum);

% A = sigma_w^2
% B = sigma_eta^2
syms A B

eqns = [ A -  E_v_sq + 2*B == 0,...
          2*B + 2*A - E_rho_sq == 0  ];
vars = [A B];

[a, b] = solve(eqns,vars);
sigma_w2 = double(a);
sigma_eta2 = double(b);

% 3) Determine optimal smoothing coefficient in exponential smoothing 

csi = sigma_w2/sigma_eta2; 
alpha = (-csi + sqrt(csi^2 + 4*csi))/2; % correct bc should be between 0,1

% 4) Perform exponential smoothing with the determined smoothing coefficient

x_hat(1) = 10;

for i = 2:n
    x_hat(i) = alpha*z(i) + (1-alpha)*x_hat(i-1);
end


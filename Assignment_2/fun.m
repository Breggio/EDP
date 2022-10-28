function [x, x_hat, z, alpha, sigma_eta2, sigma_w2, sigma] = fun(n, var_w, var_eta, incond)

rng default

a_1 = sqrt(var_w);
a_2 = sqrt(var_eta);

w = a_1.*randn(n,1);
eta = a_2.*randn(n,1);

x(1) = incond;
v_i = [];
rho_i = [];

for i = 2:n
    x(i) = x(i-1) + w(i); % generated trajectory RWM
end

for i=1:n
    z(i) = x(i) + eta(i); % Generate measurements 𝑧𝑖 of the process  𝑋𝑖
end

for i = 2:n
    v_i(i) = z(i) - z(i-1);
    if i > 2
        rho_i(i) = z(i) - z(i-2);
    end
end

E_v_sq = (1/(n-1)) *sum(v_i(1:n).^2);

E_rho_sq = (1/(n-2)) *sum(rho_i(1:n).^2);

sigma = inv([1 2; 2 2])*[E_v_sq; E_rho_sq];
sigma_w2 = sigma(1);
sigma_eta2 = sigma(2);

%3) Determine optimal smoothing coefficient in exponential smoothing 

csi = sigma_w2/sigma_eta2; 
alpha = (-csi + sqrt(csi^2 + 4*csi))/2; % correct bc should be between 0,1

% 4) Perform exponential smoothing with the determined smoothing coefficient

x_hat(1) = z(1);

for i = 2:n
    x_hat(i) = x_hat(i-1) + alpha * (z(i) - x_hat(i-1));
end

end


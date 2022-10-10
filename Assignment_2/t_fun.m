function [x_sin, z_4, x_hat_run_4] = t_fun(T,sigma_w2, sigma_eta2, a, n_4, M_4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega = 2*pi/T;
w = sqrt(sigma_w2).*randn(n_4,1);
A(1) = a;

for i = 2:n_4
    A(i) = A(i-1) + w(i);
end

x_sin = [];
for i = 1:n_4
    x_sin(i) = A(i) * sin(omega*i + 3);
end

eta = sqrt(sigma_eta2).*randn(n_4,1);

z_4 = [];
for i = 1:n_4
    z_4(i) = x_sin(i) + eta(i);
end

x_hat_run_4 = movmean(z_4, M_4);

end


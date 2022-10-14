function [z,x] = traj_meas(size_, T, sigma2_a, sigma2_n)

x = zeros(1,size_); %true data
    V = zeros(1,size_); %velocity
    z = zeros(1,size_); %measurments
    x(1) = 5; V(1) = 1;

    measurement_noise = randn*sqrt(sigma2_n);
    z(1) = x(1) + measurement_noise;
    for i=2:size_
        acceleration = randn*sqrt(sigma2_a); %randon noise of true data
        measurement_noise = randn*sqrt(sigma2_n); %random noise of measurments
        
        V(i) = V(i-1) + acceleration*T;
        x(i) = x(i-1) + V(i-1)*T + 0.5*acceleration*T^2;
        z(i) = x(i) + measurement_noise;
    end
end


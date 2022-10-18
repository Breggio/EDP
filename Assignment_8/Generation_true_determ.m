function [True_polar, True_Cartesian, Z_c, Z_p] = ...
    Generation_true_determ(N, T, sigma_D, sigma_b, InititalState)
%% Generation of true trajectory
    X   = zeros(1, N); X(1)   = InititalState(1);
    Y   = zeros(1, N); Y(1)   = InititalState(3);
    V_x = zeros(1, N); V_x(1) = InititalState(2);
    V_y = zeros(1, N); V_y(1) = InititalState(4);
    for i = 2:N
        X(i) = X(i - 1) + V_x(i - 1)*T;
        V_x(i) = V_x(i - 1);
        Y(i) = Y(i - 1) + V_y(i - 1)*T;
        V_y(i) = V_y(i - 1);
    end
    True_Cartesian = [X; V_x; Y; V_y];
    %Generation of true values of range D and aimuth b
    D = sqrt(X.^2 + Y.^2);
    b = atan(X./Y);
    True_polar = [D; b];
    %Generation of measurements
    D_m = zeros(1, N); %array of measurements of range
    b_m = zeros(1, N); %array of measurements of azimuth
    x_m = zeros(1, N); %array of pseudo-measurements of x
    y_m = zeros(1, N); %array of pseudo-measurements of y
    Z_c = zeros(2, N); %array of pseudo-measurements in Cartesian coordinates
    Z_p = zeros(2, N); %array of measurements in polar coordinates
    for i = 1:N
        D_m(i) = D(i) + randn*sigma_D;
        b_m(i) = b(i) + randn*sigma_b;

        x_m(i) = D_m(i)*sin(b_m(i));
        y_m(i) = D_m(i)*cos(b_m(i));

        Z_p(:,i) = [D_m(i); b_m(i)];
        Z_c(:,i) = [x_m(i); y_m(i)];
    end
end
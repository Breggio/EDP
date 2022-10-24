function [polar,cartesian, z_polar] = true_traj(N,T,sigma_D,sigma_b, InititalState, sigma_a, sigma_beta_add)
% Generation of true trajectory

    X   = zeros(1,N); X(1)   = InititalState(1);
    Y   = zeros(1,N); Y(1)   = InititalState(3);
    V_x = zeros(1,N); V_x(1) = InititalState(2);
    V_y = zeros(1,N); V_y(1) = InititalState(4);
    a_x = zeros(1,N - 1);
    a_y =zeros(1,N - 1);
    
    for i=2:N
        a_x(i-1) = randn * sigma_a;
        a_y(i-1) = randn * sigma_a;
        V_x(i)=V_x(i-1) + a_x(i-1)*T;
        V_y(i)=V_y(i-1) + a_y(i-1)*T;
        X(i)=X(i-1)+V_x(i-1)*T + 0.5*a_x(i-1)*T^2;
        Y(i)=Y(i-1)+V_y(i-1)*T + 0.5*a_y(i-1)*T^2;
    end
    cartesian = [X; V_x; Y; V_y];
   
    %Generation of true values of range D and aimuth b
    D = sqrt(X.^2+Y.^2);
    b = atan(X./Y);
    polar = [D;b];
    
    %Generation of measurements
    D_m=zeros(1,N); %array of measurements of range
    b_m=zeros(1,N); %array of measurements of azimuth
    
    for i=1:2:N-1
        D_m(i)=D(i)+randn*sigma_D;
        b_m(i)=b(i)+randn*sigma_b;
    end
    
    for i=4:2:N    
        D_m(i)=NaN;
        b_m(i)=b(i)+randn*sigma_beta_add;
    end
    
    z_polar = [D_m; b_m];

end
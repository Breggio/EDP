function [Z,true_data] = DataGenerator(InititalState, size_, T, sigma_a, sigma_n)
    X   = zeros(1,size_); X(1)   = InititalState(1);
    Y   = zeros(1,size_); Y(1)   = InititalState(3);
    V_x = zeros(1,size_); V_x(1) = InititalState(2);
    V_y = zeros(1,size_); V_y(1) = InititalState(4);
    a_x = zeros(1,size_ - 1);
    a_y = zeros(1,size_ - 1);
    noise_X = randn * sigma_n;
    noise_Y = randn * sigma_n;
    Z = zeros(2,size_); Z(:,1) = [X(1) + noise_X; Y(1) + noise_Y];
    for i=2:size_
        a_x(i-1) = randn * sigma_a;
        a_y(i-1) = randn * sigma_a;
        V_x(i)=V_x(i-1) + a_x(i-1)*T;
        V_y(i)=V_y(i-1) + a_y(i-1)*T;
        X(i)=X(i-1)+V_x(i-1)*T + 0.5*a_x(i-1)*T^2;
        Y(i)=Y(i-1)+V_y(i-1)*T + 0.5*a_y(i-1)*T^2;
        noise_X = randn * sigma_n;
        noise_Y = randn * sigma_n;
        Z(:,i) = [X(i)+noise_X; Y(i) + noise_Y];
    end
    true_data = [X; Y];
end


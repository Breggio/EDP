function [Z_f, X_forecast, P, K] = Kalman_Filter_2(z, m, sigma2_a, X0, P0)
    size_ = length(z);
    T = 1; sigma2_n=20^2;
    X = zeros(2, 1, size_ + 1);  %State vectors of real data
    P = zeros(2, 2, size_ + 1);  %Initial filtration error covariance matrix
    Fi = [1 T; 0 1];
    G = [0.5*T^2; T]; 
    H = [1 0]; %H

    X(:,:, 1) = X0; %Initian state vector
    P(:, :, 1) = P0; 
    Q = sigma2_a*(G*G'); %Covariance matrix of state noise
    R = sigma2_n; %Covariance matrix of measurements noise

    Z_f = zeros(1, size_);
    X_forecast = zeros(1, size_);
    K = zeros(2, 1, size_);
    for i = 2:size_ + 1
        %Prediction part
        X_pred = Fi*X(:,:,i-1);
        temp =(Fi^(m-1))*X(:,:,i-1);
        X_forecast(1,i-1) = temp(1, :);
        P_pred = ...
            Fi*P(:,:,i - 1)*Fi' + Q;
        %Filtrarion part    
        K(:,:,i-1) = P_pred*(H') * ...
            (H*P_pred*H' + R)^(-1);
        X(:,:,i) = X_pred + K(:,:,i-1)*(z(i-1) - H*X_pred);
        Z_f(i-1) = X(1,1,i);
        P(:, :, i) = (eye(2)-K(:,:,i-1)*H)*P_pred;
         
%         Use for the 15th point, where K_understimated is given  
%         %Filtrarion part    
%         X(:,:,i) = X_pred + K*(z(i-1) - H*X_pred);
%         Z_f(i-1) = X(1,1,i);
%         P(:, :, i) = (eye(2)-K*H)*P_pred;
    end
end


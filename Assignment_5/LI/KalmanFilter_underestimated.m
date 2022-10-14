function [Z_f, X_forecast, P, K] = KalmanFilter_underestimated(z, m, sigma2_a, X0, P0, K)
    N = length(z);
    T = 1; 
    sigma2_n = 20^2;
    X = zeros(2, 1, N + 1);  %State vectors of real data
    P = zeros(2, 2, N + 1);  %Initial filtration error covariance matrix
    Fi = [1 T; 0 1];
    G = [0.5*T^2; T];
    H = [1 0]; 

    X(:,:, 1) = X0; %Initian state vector
    P(:, :, 1) = P0; 
    Q = sigma2_a*(G*G'); %Covariance matrix of state noise
%     CovMatr_MeasNoise = sigma2_n; %Covariance matrix of measurements noise

    Z_f = zeros(1, N);
    X_forecast = zeros(1, N);
    for i = 2:N + 1
        %Prediction part
        X_pred = Fi*X(:,:,i-1);
        temp =(Fi^(m-1))*X(:,:,i-1);
        X_forecast(1,i-1) = temp(1, :);
        P_pred = Fi*P(:,:,i - 1)*Fi' + Q;
        %Filtrarion part    
        X(:,:,i) = X_pred + K*(z(i-1) - H*X_pred);
        Z_f(i-1) = X(1,1,i);
        P(:, :, i) = (eye(2)-K*H)*P_pred;
    end
end


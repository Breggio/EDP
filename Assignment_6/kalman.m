function [Z_f, X_f, P_pred, P, K] = kalman(z, sigma2_a, X0, P0)
%   Description:
%   CovMatr = [P; P_pred]
%   P_pred  = 2x2x(size_ + 1)
%   P = 2x2x(size_ + 1)
    
    N = length(z);
    T = 1; sigma2_n=20^2;
    Z_f = zeros(2, N + 1);                  %State vectors of real data
    P = zeros(2, 2, N + 1);                 %Filtration error covariance matrix
    P_pred = zeros(2, 2, N + 1);            %Prediction error covariance matrix
    Fi = [1 T; 0 1];
    G = [0.5*T^2; T];
    H = [1 0]; 

    Z_f(:, 1) = X0;      %Initian state vector
    P(:, :, 1) = P0; 
    Q = sigma2_a*(G*G'); %Covariance matrix of state noise
    R = sigma2_n;        %Covariance matrix of measurements noise

    X_f = zeros(1, N);
    K = zeros(2, 1, N);
    for i = 2:N + 1
        %Prediction part
        X_pred = Fi*Z_f(:,i-1);
        P_pred(:,:,i-1) = Fi * P(:,:,i - 1) * Fi' + Q;
        %Filtrarion part    
        K(:,:,i-1) = P_pred(:,:,i-1)*(H') * (H*P_pred(:,:,i-1)*H' + R)^(-1);
        Z_f(:,i) = X_pred + K(:,:,i-1)*(z(i-1) - H*X_pred);
        P(:, :, i) = (eye(2)-K(:,:,i-1)*H)*P_pred(:,:,i-1);
    end
end


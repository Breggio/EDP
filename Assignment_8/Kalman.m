function [Z_filtered, P, P_pred, range_fe, azimuth_fe, CondMatr, K] ...
    = Kalman(Z_cart, Z_polar, T, sigma_D, sigma_beta)
%   Description:
    N = length(Z_cart);
    Z_filtered = zeros(4, N + 1); %Filtered data
    Z_forecast = zeros(4, N);     %Forecast data
    P = zeros(4, 4, N + 1);       %Filtration error covariance matrix
    P_pred = zeros(4, 4, N + 1);  %Prediction error covariance matrix
    Fi = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];
    H = [1 0 0 0; 0 0 1 0];
    Z_filtered(:, 1) = [40000; -20; 40000; -20]; %Initian state vector
    P(:, :, 1) = [10^10 0 0 0; 0 10^10 0 0; 0 0 10^10 0; 0 0 0 10^10];
    R = zeros(2,2);
    K = zeros(4, 2, N);
    range_fe = zeros(2,N);   %array of filtered and extrapolated range
    azimuth_fe = zeros(2,N); %array of filtered and extrapolated azimuth
    CondMatr = zeros(1,N);
    
    for i = 2:N + 1
        R(1,1) = sin(Z_polar(2,i - 1))^2*sigma_D^2 + (Z_polar(1,i - 1))^2*cos(Z_polar(2,i - 1))^2*sigma_beta^2;
        R(2,2) = cos(Z_polar(2,i - 1))^2*sigma_D^2 + (Z_polar(1,i - 1))^2*sin(Z_polar(2,i - 1))^2*sigma_beta^2;
        R(1,2) = sin(Z_polar(2,i - 1))*cos(Z_polar(2,i - 1))*(sigma_D^2 - (Z_polar(1,i - 1))^2*sigma_beta^2);
        R(2,1) = sin(Z_polar(2,i - 1))*cos(Z_polar(2,i - 1))*(sigma_D^2 - (Z_polar(1,i - 1))^2*sigma_beta^2);
        %Prediction part
        Z_forecast(:,i - 1) =  Fi * Z_filtered(:,i - 1);
        P_pred(:,:,i-1) = Fi * P(:,:,i - 1) * Fi';
        %Filtrarion part    
        K(:,:,i - 1) = P_pred(:,:,i - 1)*(H') * (H*P_pred(:,:,i - 1)*H' + R)^(-1);
        P(:, :, i) = (eye(4) - K(:,:,i - 1)*H)*P_pred(:,:,i - 1);
        Z_filtered(:,i) = Z_forecast(:,i - 1) + K(:,:,i - 1)*(Z_cart(:,i - 1) - H*Z_forecast(:,i - 1));
        
        range_fe(1,i - 1) = sqrt( Z_filtered(1,i)^2 + Z_filtered(3,i)^2 );
        range_fe(2,i - 1) = sqrt( Z_forecast(1,i - 1)^2 + Z_forecast(3,i - 1)^2);
        
        azimuth_fe(1,i - 1) = atan(Z_filtered(1,i)/Z_filtered(3,i));
        azimuth_fe(2,i - 1) = atan(Z_forecast(1,i - 1)/Z_forecast(3,i - 1));
        
        if (sigma_D^2) > (Z_polar(1,i - 1)^2*sigma_beta^2)
            CondMatr(i - 1) = (sigma_D^2)/(Z_polar(1,i - 1)^2*sigma_beta^2);
        else
            CondMatr(i - 1) = (Z_polar(1,i - 1)^2*sigma_beta^2)/(sigma_D^2);
        end
    end
end


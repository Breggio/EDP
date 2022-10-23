function [Z_filtered, FiltrErr_CovMatr, ForecErr_CovMatr] = Kalman_v2(Z, T, sigma_a, sigma_n, P0, Z0)
    size_ = length(Z);
    TransMatr = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];
    ObservMatr = [1 0 0 0; 0 0 1 0];
    InputMatr = [0.5*T^2 0; T 0; 0 0.5*T^2; 0 T];
    sigma_d = 200^2;
    sigma_beta = 0.01^2;
    CovMatr_StateNoise = (InputMatr*InputMatr.')*sigma_a^2;
    CovMatr_MeasNoise = [40000 0; 0 0.01^2];

    
    Z_filtered = zeros(4, size_ + 1);           %Filtered data
    Z_forecast = zeros(4, size_);               %Forecast data
    FiltrErr_CovMatr = zeros(4, 4, size_ + 1);  %Filtration error covariance matrix
    ForecErr_CovMatr = zeros(4, 4, size_);      %Prediction error covariance matrix
    K = zeros(4, 2, size_);
    
    Z_filtered(:, 1) = Z0;
    FiltrErr_CovMatr(:, :, 1) = P0;
    for i = 2:size_+1
        %Prediction part
        Z_forecast(:,i-1) =  TransMatr*Z_filtered(:,i-1);
        ForecErr_CovMatr(:,:,i-1) = ...
            TransMatr*FiltrErr_CovMatr(:,:,i - 1)*TransMatr.' + CovMatr_StateNoise;
        %Filtrarion part    
        K(:,:,i-1) = ForecErr_CovMatr(:,:,i-1)*(ObservMatr.') * ...
            (ObservMatr*ForecErr_CovMatr(:,:,i-1)*ObservMatr.' + CovMatr_MeasNoise)^(-1);
        FiltrErr_CovMatr(:, :, i) = ...
            (eye(4)-K(:,:,i-1)*ObservMatr)*ForecErr_CovMatr(:,:,i-1);
        Z_filtered(:,i) = Z_forecast(:,i-1) + ...
            K(:,:,i-1)*(Z(:,i-1) - ObservMatr*Z_forecast(:,i-1));
    end
    Z_filtered = Z_filtered(:,2:size_ + 1);
    FiltrErr_CovMatr = FiltrErr_CovMatr(:,:,2:size_ + 1);
end

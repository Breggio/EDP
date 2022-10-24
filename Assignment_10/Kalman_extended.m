function [Z_filtered, Z_forecast, range_fe, azimuth_fe] ...
    = Kalman_extended(Z_polar, T, sigma_D,sigma_b, sigma_a, sigma_b_add)

%   Description:
    size_ = length(Z_polar);
    Z_filtered = zeros(4, size_);           %Filtered data
    Z_forecast = zeros(4, size_);           %Forecast data [x; Vx; y; Vy]
    FiltrErr_CovMatr = zeros(4, 4, size_);  %Filtration error covariance matrix
    PredErr_CovMatr  = zeros(4, 4, size_); %Prediction error covariance matrix
    K = zeros(4, 2, size_);
    
    Z_filtered(:, 3) = ...
        [Z_polar(1,3)*sin(Z_polar(2,3));...
        (Z_polar(1,3)*sin(Z_polar(2,3)) - Z_polar(1,1)*sin(Z_polar(2,1)))/(2*T);...
        Z_polar(1,3)*cos(Z_polar(2,3));...
        (Z_polar(1,3)*cos(Z_polar(2,3)) - Z_polar(1,1)*cos(Z_polar(2,1)))/(2*T)]; %Initian state vector
    FiltrErr_CovMatr(:, :, 3) = [10^4 0 0 0; 0 10^4 0 0;...
        0 0 10^4 0; 0 0 0 10^4];
    
    TransMatr = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1]; %Ô
    InputMatr = [0.5*T^2 0; T 0; 0 0.5*T^2; 0 T]; %G
    CovMatr_StateNoise   = (InputMatr*InputMatr')*sigma_a^2; %Covariance matrix of state noise

    range_fe = zeros(2,size_); %array of filtered and extrapolated range
    azimuth_fe = zeros(2,size_); %array of filtered and extrapolated azimuth
    for i = 4:1:size_
        %Forecasting
        Z_forecast(:,i-1) =  TransMatr*Z_filtered(:,i-1);
        PredErr_CovMatr(:,:,i-1) = ...
            TransMatr*FiltrErr_CovMatr(:,:,i-1)*TransMatr' + CovMatr_StateNoise;
        h = [sqrt(Z_forecast(1,i-1)^2 + Z_forecast(3,i-1)^2 );...
             atan(Z_forecast(1,i-1)/Z_forecast(3,i-1)) ];
        dhdx = zeros(2, 4);
        dhdx(1,1) =  Z_forecast(1,i-1)/sqrt( (Z_forecast(1,i-1))^2 + (Z_forecast(3,i-1))^2 );
        dhdx(1,3) =  Z_forecast(3,i-1)/sqrt( (Z_forecast(1,i-1))^2 + (Z_forecast(3,i-1))^2 );
        dhdx(2,1) =  Z_forecast(3,i-1)/(Z_forecast(1,i-1)^2 + Z_forecast(3,i-1)^2 );
        dhdx(2,3) = -Z_forecast(1,i-1)/(Z_forecast(1,i-1)^2 + Z_forecast(3,i-1)^2 );
        
        if mod(i,2) == 0
            CovMatr_MeasureNoise = [sigma_D^2 0; 0 sigma_b_add^2];
        else
            CovMatr_MeasureNoise = [sigma_D^2 0; 0 sigma_b^2];
        end
        
        %Filtrarion part 
        K(:,:,i) = PredErr_CovMatr(:,:,i-1)*(dhdx') * ...
            (dhdx*PredErr_CovMatr(:,:,i-1)*dhdx' + CovMatr_MeasureNoise)^(-1);
        
        FiltrErr_CovMatr(:, :, i) = (eye(4)-K(:,:,i)*dhdx)*PredErr_CovMatr(:,:,i-1);
        Z_temp = Z_polar(:,i);
        if isnan(Z_temp(1))
            Z_temp(1) = h(1);
        end
        Z_filtered(:,i) = Z_forecast(:,i - 1) + ...
            K(:,:,i)*(Z_temp - h);

        range_fe(1,i)=sqrt( Z_filtered(1,i)^2+Z_filtered(3,i)^2 );
        range_fe(2,i)=sqrt( Z_forecast(1,i-1)^2+Z_forecast(3,i-1)^2);
        azimuth_fe(1,i)=atan(Z_filtered(1,i)/Z_filtered(3,i));
        azimuth_fe(2,i)=atan(Z_forecast(1,i-1)/Z_forecast(3,i-1));
    end
    azimuth_fe(1,1:3) = azimuth_fe(1,4);
end

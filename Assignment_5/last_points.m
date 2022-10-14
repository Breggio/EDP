clc
clear
%% Point 11
size_ = 200; T = 1;
X0 = [2; 0]; P0 = [10000 0; 0 10000];
sigma2_n = 20^2; sigma2_a = 0.2^2;
m = 7; M = 500;
ErrFiltered = zeros(M, 1, size_);

for i = 1:M
    [z, x] = DataGenerator(size_, T, sigma2_a, sigma2_n);
    [Z_filtered, X_forecast, P_matr,] =...
        KalmanFilter(z, m, sigma2_a, X0, P0);
    ErrFiltered(i,1,:) = (x - Z_filtered).^2;
end
Final_ErrFiltered = sqrt( 1/(M-1) * sum(ErrFiltered) ); %True estimation error
sqrt_CovMatr_ErrFiltr = sqrt(P_matr(1,1,2:size_ + 1));

figure(11)
plot(3:size_, Final_ErrFiltered(1,3:200), 'b', 3:size_, sqrt_CovMatr_ErrFiltr(1, 3:200), 'r')
title('Point 11. Compare calculation errors of estimation and true estimation errors.')
legend('True estimation error', 'Standard deviation of estimation error')
xlabel('Observation interval')
ylabel('True estimation error | SD of estimation error')

%% Pont 12
sigma2_a = 0;

for i = 1:M
    [z, x] = DataGenerator(size_, T, sigma2_a, sigma2_n);
    [Z_filtered, X_forecast, P_matr, K_arr] = ...
        KalmanFilter(z, m, sigma2_a, X0, P0);
    ErrFiltered(i,1,:) = (x - Z_filtered).^2;    
end
Final_ErrFiltered = sqrt( 1/(M-1) * sum(ErrFiltered) );
sqrt_CovMatr_ErrFiltr = sqrt(P_matr(1,1,2:size_ + 1));

figure(12)
ax1 = subplot(2,2,1:2);
plot(1:size_, x, 'b', 1:size_, z, 'r.', 1:size_, Z_filtered, 'g')
legend('True','Measurements','Filtered process')
title(ax1, 'Process example')
xlabel('Observation interval')
ylabel('Data')

ax2 = subplot(2,2,3);
plot(1:size_ , K_arr(1,:))
legend('Gain')
title(ax2, 'Filter gain')
xlabel('Observation interval')
ylabel('Gain value')

ax3 = subplot(2,2,4);
plot(3:size_, Final_ErrFiltered(1,3:size_),...
    3:size_, sqrt_CovMatr_ErrFiltr(1, 3:200))
legend('True estimation errors', 'Calculation errors of estimation')
title(ax3, 'Errors comparing')
xlabel('Observation interval')
ylabel('True error and SD of estimation error')

%% Pont 13
ErrForecast = zeros(M, 1, size_);

for i = 1:M
    [z, x] = DataGenerator(size_, T, 0.2^2, sigma2_n); %sigma2_a = 0.2^2;
    [Z_filtered, X_forecast, P_matr, K_arr] = ...
        KalmanFilter(z, m, 0, X0, P0); %sigma2_a = 0;
    ErrFiltered(i,1,:) = (x - Z_filtered).^2;
    ErrForecast(i,1,:) = (x - X_forecast).^2;
end
Final_ErrFiltered = sqrt( 1/(M-1) * sum(ErrFiltered) );
Final_ErrForecast = sqrt( 1/(M-1) * sum(ErrForecast) );

figure(13)
ax1 = subplot(2,2,1:2);
plot(1:size_, x, 'b', 1:size_, z, 'r.', 1:size_, Z_filtered, 'g')
legend('True','Measurements','Filtered process')
title(ax1, 'Process example')
xlabel('Observation interval')
ylabel('Data')

ax2 = subplot(2,2,3);
plot(3:size_, Final_ErrFiltered(1,3:size_),...
     3:size_, Final_ErrForecast(1,3:size_))
legend('Filtered estimate error', 'Forecast estimate error')
title(ax2, 'Estimate errors comparing')
xlabel('Observation interval')
ylabel('Filtered and Forecast estimates')

ax3 = subplot(2,2,4);
plot(3:size_, Final_ErrFiltered(1,3:size_),...
    3:size_, sqrt_CovMatr_ErrFiltr(1, 3:200))
legend('True estimation errors', 'Calculation errors of estimation')
title(ax3, 'Errors comparing')
xlabel('Observation interval')
ylabel('True error and SD of estimation error')

%% Point 14
sigma2_a = 1;
[z, x] = DataGenerator(size_, T, sigma2_a, sigma2_n);
[Z_filtered, X_forecast, P_matr, K_arr] = ...
        KalmanFilter(z, m, sigma2_a, X0, P0);

figure(14)
ax1 = subplot(2,2,1);
plot(1:size_, x, 'b', 1:size_, z, 'r.', 1:size_, Z_filtered, 'g')
legend('True','Measurements','Filtered process')
title(ax1, 'Process with \sigma_a = 1')
xlabel('Observation interval')
ylabel('Data')

ax3 = subplot(2,2,3);
plot(1:size_ , K_arr(1,:))
legend('Gain')
title(ax3, 'Filter gain with \sigma_a = 1')
xlabel('Observation interval')
ylabel('Gain value')
grid on

sigma2_a = 0.2^2;
[z, x] = DataGenerator(size_, T, sigma2_a, sigma2_n); 
[Z_filtered, X_forecast, P_matr, K_arr] = ...
        KalmanFilter(z, m, sigma2_a, X0, P0);

ax2 = subplot(2,2,2);
plot(1:size_, x, 'b', 1:size_, z, 'r.', 1:size_, Z_filtered, 'g')
legend('True','Measurements','Filtered process')
title(ax2, 'Process with \sigma_a = 0.2')
xlabel('Observation interval')
ylabel('Data')

ax4 = subplot(2,2,4);
plot(1:size_ , K_arr(1,:))
legend('Gain ')
title(ax4, 'Filter gain with \sigma_a = 0.2')
xlabel('Observation interval')
ylabel('Gain value')
grid on

%% Point 15
X0 = [100; 5];

for i = 1:M
    [z, x] = DataGenerator(size_, T, sigma2_a, sigma2_n); %sigma2_a = 0.2^2;
    [Z_filtered, X_forecast, P_matr, K_arr] = ...
        KalmanFilter(z, m, sigma2_a, X0, P0); %sigma2_a = 0;
    ErrFiltered(i,1,:) = (x - Z_filtered).^2;
end
Final_ErrFiltered = sqrt( 1/(M-1) * sum(ErrFiltered) );

figure(15)
ax1 = subplot(2,2,1);
plot(1:size_, x, 'b', 1:size_, z, 'r.', 1:size_, Z_filtered, 'g')
legend('True','Measurements','Filtered process')
title(ax1, 'Process example with optimal gain')
xlabel('Observation interval')
ylabel('Data')

ax3 = subplot(2,2,3:4);
plot(3:size_, Final_ErrFiltered(1,3:size_))
title(ax3, 'Estimate errors comparing')
xlabel('Observation interval')
ylabel('Filtered estimates')

K_underestimated = K_arr(1,size_)/1;
for i = 1:M
    [z, x] = DataGenerator(size_, T, sigma2_a, sigma2_n); %sigma2_a = 0.2^2;
    [Z_filtered, X_forecast, P_matr, ] = ...
        KalmanFilter_underestimated(z, m, sigma2_a, X0, P0, K_underestimated); %sigma2_a = 0;
    ErrFiltered(i,1,:) = (x - Z_filtered).^2;
end
Final_ErrFiltered = sqrt( 1/(M-1) * sum(ErrFiltered) );

ax2 = subplot(2,2,2);
plot(1:size_, x, 'b', 1:size_, z, 'r.', 1:size_, Z_filtered, 'g')
legend('True','Measurements','Filtered process underestimated gain')
title(ax2, 'Process example with underestimated gain')
xlabel('Observation interval')
ylabel('Data')

subplot(2,2,3:4); hold on
plot(3:size_, Final_ErrFiltered(1,3:size_))
legend('Filtered estimate error with optimal filter gain',...
    'Forecast estimate error underestimated filter gain')
grid on
hold off
text(150, 70 , "K optimal = " + num2str(K_arr(1,size_)))
text(150, 60 , "K underestimated = " + num2str(K_underestimated))
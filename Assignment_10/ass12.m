%% Joint Assimilation of Navigation Datan Coming from Different Sources

% Written by Irina Yareshko and Luca Breggion, Skoltech 2022

close all 
clear
clc

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 

T=2; 
N=500;

InititalState = [1000;100;1000;100]; %[x0; Vx0; y0; Vy0]

sigma_a = 0.3;          % variance of acceleration noise
sigma_D = 50;           % variance of range noise of measurements
sigma_b = 0.004;        % variance of azimuth noise of measurements
sigma_b_add = 0.001;    %

% Creating of arrays for running Kalman filter M times

M=500;                            %number of runs
err_range_filt  = zeros(M,N); %Filtraion error of range
err_range_pred  = zeros(M,N); %Predicted error of range
err_azimut_filt = zeros(M,N); %Filtraion error of azimut
err_azimut_pred = zeros(M,N); %Predicted error of azimut

for i=1:M
    %Generation of deterministic trajectory and its measurements
    [polar,cart,z_p] = ...
        true_traj(N,T,sigma_D,sigma_b,InititalState, sigma_a, sigma_b_add);

    %Applying or Kalman filter for measurements
    [Z_filtered,Z_forecast, range_fe, azimuth_fe] = Kalman_extended(z_p,...
        T,sigma_D,sigma_b, sigma_a, sigma_b_add);
    
    err_range_filt(i,:) = (polar(1,:) - range_fe(1,:)).^2;
    err_range_pred(i,:) = (polar(1,:) - range_fe(2,:)).^2;
    err_azimut_filt(i,:)= (polar(2,:) - azimuth_fe(1,:)).^2;
    err_azimut_pred(i,:)= (polar(2,:) - azimuth_fe(2,:)).^2;
end

err_range_filt =sqrt(1/(M-1)*sum(err_range_filt));
err_range_pred =sqrt(1/(M-1)*sum(err_range_pred));
err_azimut_filt = sqrt(1/(M-1)*sum(err_azimut_filt));
err_azimut_pred = sqrt(1/(M-1)*sum(err_azimut_pred));

%% Point 6

% True trajectory, measurements and filtered&extrapolated
figure(1)
polarplot(polar(2,:),polar(1,:), 'm',...
    z_p(2,:),z_p(1,:),'c.', ...
    azimuth_fe(1,:), range_fe(1,:), 'k')
legend('True trajectory', 'Measurements', 'Filtered', 'FontSize', 30)
grid on; grid minor

%% Point 7

figure(2)
plot(4:N,err_range_filt(4:N),'m',...
     4:N,err_range_pred(4:N),'c');
legend('Filtration error','Extrapolation error');
%title('Errors of range');
xlabel('Step')
ylabel('Errors')
grid on; grid minor

%% Point 7

figure(3)
plot(4:N,err_azimut_filt(4:N),'m', ...
     4:N,err_azimut_pred(4:N),'c');
legend('Filtration error','Extrapolation error');
%title('Errors of azimiuth');
xlabel('Step')
ylabel('Errors')
grid on; grid minor

%% Point 7

figure(4)
plot(cart(1,:),cart(3,:), 'm',...
    Z_filtered(1,:),Z_filtered(3,:),'c')
%title(True and filtered trajectory')
legend('True', 'Filtered')
grid on; grid minor
xlabel('Step')
ylabel('Data')

%% Point 8
% A = NaN;
% b = z_p(1,:);
% b(isnan(A))=0;

z_pp = z_p(1,1:2:N);

figure(7)
plot(1:N, polar(1,:),'m', 1:2:N, z_pp, 'c',1:N ,range_fe(1,:), 'k')
%title('Range')
legend('True', 'Measurements', 'Filtered', 'location', 'best','FontSize', 25)
grid on; grid minor
xlabel('Step')
ylabel('Data')

%% Point 8

figure(8)
plot(1:N, polar(2,1:N),'m', 4:N, z_p(2,4:N), 'c', 1:N, azimuth_fe(1,:), 'k')
%title('Azimuth')
legend('True', 'Measurements', 'Filtered')
grid on; grid minor
xlabel('Step')
ylabel('Data')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                FUNCTION                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


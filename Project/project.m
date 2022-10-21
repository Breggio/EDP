%% Estimation of a site where motion of a moving vehicle started using radar data
% Coordinate transformation of measurements 

% Written by Irina Yareshko and Luca Breggion, Skoltech 2022

close all 
clear
clc

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 2)

z1 = importdata('z1.txt'); % z1
z2 = importdata('z2.txt'); % z2
z3 = importdata('z3.txt'); % z3

% First column - measurements of the range 𝐷  in meters 
% Second column - measurements of the azimuth 𝛽  in radians 

sigma_a = 0.01^2;
sigma_d = 200^2;
sigma_beta = 0.01^2;
N = length(z1);
T = 1;
M = 500;

P0 = 10^4 * eye(4);

in_state = zeros(4,N);

X = zeros(1,N);
Y = zeros(1,N);
V_x = zeros(1,N);
V_y = zeros(1,N);

X(1) = in_state(1);
Y(1) = in_state(3);
V_x(1) = in_state(2);
V_y(1) = in_state(4);

cnd_n = zeros(1,N);
cnd_matr = zeros(1,N);

d = [z1(:,1), z2(:,1), z3(:,1)];
beta = [z1(:,2), z2(:,2), z3(:,2)];

for i = 1:N
    a=randn;
    d_m_1(i) = d(i,1) + a*sigma_d;
    b_m_1(i) = beta(i,1) + a*sigma_beta;
    d_m_2(i) = d(i,2) +a*sigma_d;
    b_m_2(i) = beta(i,2) + a*sigma_beta;
    d_m_3(i) = d(i,3) + a*sigma_d;
    b_m_3(i) = beta(i,3) + a*sigma_beta;

    x_m_1(i) = d_m_1(i)*sin(b_m_1(i));
    y_m_1(i) = d_m_1(i)*cos(b_m_1(i));
    x_m_2(i) = d_m_2(i)*sin(b_m_2(i));
    y_m_2(i) = d_m_2(i)*cos(b_m_2(i));
    x_m_3(i) = d_m_3(i)*sin(b_m_3(i));
    y_m_3(i) = d_m_3(i)*cos(b_m_3(i));

    Z_p_1(:,i) = [d_m_1(i); b_m_1(i)];
    Z_p_2(:,i) = [d_m_2(i); b_m_2(i)];
    Z_p_3(:,i) = [d_m_3(i); b_m_3(i)];

    Z_c_1(:,i) = [x_m_1(i); y_m_1(i)];
    Z_c_2(:,i) = [x_m_2(i); y_m_2(i)];
    Z_c_3(:,i) = [x_m_3(i); y_m_3(i)];

end

%%

clear z1
clear z2
clear z3

sigma_n = 20;
z1 = Z_c_1;
z2 = Z_c_2;
z3 = Z_c_3;

sigma_a = 0.01^2;
sigma_n = 20;

size_ = N;

Z0 = [0;0;0;0]; 
X0 = [2; 0; 2; 0]; 

SqrError_RAW_f = zeros(M, 1, size_ - 2); %Squared error of coordinate
SqrError_RAW_s = zeros(M, 1, size_ - 2); %Squared error of coordinate
SqrError_OPT_f = zeros(M, 1, size_ - 2); %Squared error of velocity
SqrError_OPT_s = zeros(M, 1, size_ - 2); %Squared error of velocity

for i = 1:M
    [Z_XY, True_XY] = DataGenerator(X0, size_, T, sigma_a, sigma_n);
    P0 = 10^4 * eye(4);
    
    [z_raw_f, FiltrErr_CovMatr, ForecErr_CovMatr] = Kalman_v2(Z_XY, T, sigma_a, sigma_n, P0, Z0);
    [z_raw_s, ] = Kalman_BackSmooth(z_raw_f, ForecErr_CovMatr, FiltrErr_CovMatr, T);
    
    P0 = FiltrErr_CovMatr(:,:,size_);
    [z_opt_f, FiltrErr_CovMatr, ForecErr_CovMatr] = Kalman_v2(Z_XY, T, sigma_a, sigma_n, P0, Z0);
    [z_opt_s, ] = Kalman_BackSmooth(z_opt_f, ForecErr_CovMatr, FiltrErr_CovMatr, T);

    SqrError_RAW_f(i,:) = ( True_XY(1,3:size_) - z_raw_f(1,3:size_) ).^2;
    SqrError_RAW_s(i,:) = ( True_XY(1,3:size_) - z_raw_s(1,3:size_) ).^2;
    SqrError_OPT_f(i,:) = ( True_XY(1,3:size_) - z_opt_f(1,3:size_) ).^2;
    SqrError_OPT_s(i,:) = ( True_XY(1,3:size_) - z_opt_s(1,3:size_) ).^2;
end

Final_ErrSmoothed_RAW_f = sqrt( 1/(M-1) * sum(SqrError_RAW_f) ); %true estimation error of cordinate
Final_ErrSmoothed_RAW_s = sqrt( 1/(M-1) * sum(SqrError_RAW_s) ); %true estimation error of cordinate
Final_ErrSmoothed_OPT_f = sqrt( 1/(M-1) * sum(SqrError_OPT_f) ); %true estimation error of velocity
Final_ErrSmoothed_OPT_s = sqrt( 1/(M-1) * sum(SqrError_OPT_s) ); %true estimation error of velocity

% figure(1)
% subplot(1,2,1)
% plot(True_XY(1,:), True_XY(2,:), 'black',...
%      Z_XY(1,:), Z_XY(2,:), 'blue.',...
%      True_XY(1,1), True_XY(2,1), 'redo');
% legend('True trajectory', 'GPS Data', 'Start point')
% xlabel('X'), ylabel('Y')
% title('Example of true trajectory and GPS data distribution')
% text(40, -30, join(['\sigma_a = ', num2str(0.01)]))
% text(40, -25, join(['\sigma_n = ', num2str(20.0)]))
% grid on
% subplot(1,2,2)
% plot(z_raw_s(1,:), z_raw_s(3,:), 'b',...
%      z_opt_s(1,:), z_opt_s(3,:), 'r',...
%      True_XY(1,:), True_XY(2,:), 'black');
% legend('z raw smoothed', 'z opt smoothed', 'true trajectory')
% xlabel('X'), ylabel('Y')
% title('Example of correction')
% grid on

figure(2)
plot(3:size_, Final_ErrSmoothed_RAW_f(1,:), 'black',...
     3:size_, Final_ErrSmoothed_RAW_s(1,:), 'red',...
     3:size_, Final_ErrSmoothed_OPT_f(1,:), 'blue',...
     3:size_, Final_ErrSmoothed_OPT_s(1,:), 'green')
legend('Non-Tuned Kalman',...
       'Non-Tuned Kalman + Backward Smoothing',...
       'Tuned Kalman',...
       'Tuned Kalman + Backward Smoothing');
title('Final True Error');
xlabel('Observe interval');
ylabel('Error')
grid on

%% 

P0 = 10^4 * eye(4);
size_ = length(z1); 
sigma_a = 0.01; 
sigma_n = 15^2; 
T = 1;

[z1_filtered_raw, FiltrErr_CovMatr, ForecErr_CovMatr] = Kalman_v2(z1, T, sigma_a, sigma_n, P0, Z0);
[z1_smoothed_raw, ] = Kalman_BackSmooth(z1_filtered_raw, ForecErr_CovMatr, FiltrErr_CovMatr, T);

[z2_filtered_raw, FiltrErr_CovMatr, ForecErr_CovMatr] = Kalman_v2(z2, T, sigma_a, sigma_n, P0, Z0);
[z2_smoothed_raw, ] = Kalman_BackSmooth(z2_filtered_raw, ForecErr_CovMatr, FiltrErr_CovMatr, T);

[z3_filtered_raw, FiltrErr_CovMatr, ForecErr_CovMatr] = Kalman_v2(z3, T, sigma_a, sigma_n, P0, Z0);
[z3_smoothed_raw, ] = Kalman_BackSmooth(z3_filtered_raw, ForecErr_CovMatr, FiltrErr_CovMatr, T);

P0 = FiltrErr_CovMatr(:,:,499);
[z1_filtered, FiltrErr_CovMatr, ForecErr_CovMatr] = Kalman_v2(z1, T, sigma_a, sigma_n, P0, Z0);
[z1_smoothed, ] = Kalman_BackSmooth(z1_filtered, ForecErr_CovMatr, FiltrErr_CovMatr, T);

[z2_filtered, FiltrErr_CovMatr, ForecErr_CovMatr] = Kalman_v2(z2, T, sigma_a, sigma_n, P0, Z0);
[z2_smoothed, ] = Kalman_BackSmooth(z2_filtered, ForecErr_CovMatr, FiltrErr_CovMatr, T);

[z3_filtered, FiltrErr_CovMatr, ForecErr_CovMatr] = Kalman_v2(z3, T, sigma_a, sigma_n, P0, Z0);
[z3_smoothed, ] = Kalman_BackSmooth(z3_filtered, ForecErr_CovMatr, FiltrErr_CovMatr, T);

temp = FiltrErr_CovMatr(1,1,1:200);

figure()
plot(z1(1,:), z1(2,:), 'm*',z2(1,:), z2(2,:),'c*', z3(1,:), z3(2,:),'k*')
grid on; grid minor
legend('Vehicle 1', 'Vehicle 2', 'Vehicle 3')
title('Measurements')
xlabel('X [m]')
xlabel('Y [m]')
% xlim([0 z1(1,end)])
% ylim([0 -z1(1,end)])

%% Our plot

figure()
P = plot(z1(1,:), z1(2,:), '.',z2(1,:), z2(2,:),'.', z3(1,:), z3(2,:),'.',...
    z1_filtered_raw(1,:), z1_filtered_raw(3,:),'m', ...
    z2_filtered_raw(1,:), z2_filtered_raw(3,:),'c', ...
    z3_filtered_raw(1,:), z3_filtered_raw(3,:),'k');
P(1).LineWidth = 0.00001;
P(2).LineWidth = 0.00001;
P(3).LineWidth = 0.00001;
P(4).LineWidth = 1.5;
P(5).LineWidth = 1.5;
P(6).LineWidth = 1.5;
grid on; grid minor
title('Filtered trajectory with the initial P')
legend('Measurements Vehicle 1', 'Measurements Vehicle 2', 'Measurements Vehicle 3',...
    'Vehicle 1', 'Vehicle 2', 'Vehicle 3')
xlabel('X [m]')
ylabel('Y [m]')

figure()
P = plot(z1(1,:), z1(2,:), '.',z2(1,:), z2(2,:),'.', z3(1,:), z3(2,:),'.',...
    z1_filtered(1,:), z1_filtered(3,:),'m', ...
    z2_filtered(1,:), z2_filtered(3,:),'c', ...
    z3_filtered(1,:), z3_filtered(3,:),'k');
P(1).LineWidth = 0.00001;
P(2).LineWidth = 0.00001;
P(3).LineWidth = 0.00001;
P(4).LineWidth = 2;
P(5).LineWidth = 2;
P(6).LineWidth = 2;
grid on; grid minor
title('Filtered trajectory with the optimal P')
legend('Measurements Vehicle 1', 'Measurements Vehicle 2', 'Measurements Vehicle 3',...
    'Vehicle 1', 'Vehicle 2', 'Vehicle 3')
xlabel('X [m]')
ylabel('Y [m]')

figure()
P = plot(z1(1,:), z1(2,:), '.',z2(1,:), z2(2,:),'.', z3(1,:), z3(2,:),'.',...
    z1_smoothed(1,:), z1_smoothed(3,:),'m', ...
    z2_smoothed(1,:), z2_smoothed(3,:),'c', ...
    z3_smoothed(1,:), z3_smoothed(3,:),'k', z1_smoothed(1,1), z1_smoothed(3,1), ...
    z2_smoothed(1,1), z2_smoothed(3,1), z3_smoothed(1,1), z3_smoothed(3,1), ...
    z1_smoothed(1,499), z1_smoothed(3,499), ...
    z2_smoothed(1,499), z2_smoothed(3,499), z3_smoothed(1,499), z3_smoothed(3,499));
P(1).LineWidth = 0.00001;
P(2).LineWidth = 0.00001;
P(3).LineWidth = 0.00001;
P(4).LineWidth = 2;
P(5).LineWidth = 2;
P(6).LineWidth = 2;
grid on; grid minor
legend('Vehicle 1', 'Vehicle 2', 'Vehicle 3')
title('Smoothed Trajectory')
xlabel('X [m]')
ylabel('Y [m]')

%%

figure(3)
temp2 = FiltrErr_CovMatr(1,1,1:200);
plot(1:200, temp(1,:), 1:200, temp2(1,:))
legend('1st run. Non-tuned initial P_0', '2nd run. Tuned initial P_0')
%title('Comparing filter error covarience matrixes')
xlabel('Observation interval')
ylabel('Error SD')
grid on; grid minor

figure(4)
plot(z1(1,:), z1(2,:), 'blue*',...
     z2(1,:), z2(2,:), 'red*',...
     z3(1,:), z3(2,:), 'green*')
legend('Vehicle1','Vehicle2','Vehicle3')
title('Measurements')
xlabel('X'), ylabel('Y')
grid on; grid minor

figure(5)
plot(z1_smoothed_raw(1,:), z1_smoothed_raw(3,:), 'blue',...
     z2_smoothed_raw(1,:), z2_smoothed_raw(3,:), 'red',...
     z3_smoothed_raw(1,:), z3_smoothed_raw(3,:), 'green')
legend('Vehicle1','Vehicle2','Vehicle3')
title('Measurements')
xlabel('X'), ylabel('Y')
grid on; grid minor
% xlim([-100 500])
% ylim([-100 600])

figure(6)
plot(z1_smoothed(1,:), z1_smoothed(3,:), 'red',...
     z2_smoothed(1,:), z2_smoothed(3,:), 'blue',...
     z3_smoothed(1,:), z3_smoothed(3,:), 'green',...
     z1_smoothed(1,1), z1_smoothed(3,1), 'redo',...
     z2_smoothed(1,1), z2_smoothed(3,1), 'blueo',...
     z3_smoothed(1,1), z3_smoothed(3,1), 'greeno')
legend('Vehicle 1','Vehicle 2','Vehicle 3')
title('Corrected trajectory')
xlabel('x'), ylabel('y')
grid on
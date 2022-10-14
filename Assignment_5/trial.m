close all; clear; clc;

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 1-2)
%Creating of arrays
x=zeros(1,200); %true data
V=zeros(1,200); %velocity
Z=zeros(1,200); %measurments

sigma2_n=20^2;
n=randn*sqrt(sigma2_n); %random noise of measurments
%Initial data
N=200;
T=1;
x(1)=5;
V(1)=1;
Z(1)=x(1)+n; %the first measurment

%Generation of data
for i=2:N
    sigma2_a=0.2^2;
    a=randn*sqrt(sigma2_a); %normally distributed random acceleration
    n=randn*sqrt(sigma2_n); %random noise of measurments
    V(i)=V(i-1)+a*T;
    x(i)=x(i-1)+V(i-1)*T+a*T^2/2;
    Z(i)=x(i)+n;
end

%% 4) 

%Presenting the system at state space
%State matrixes
Fi=[1 T; 0 1]; %transition matrix
G=[T/2;T]; %input matrix
H=[1 0]; %observation matrix

%Kalman filter development
%Creating of arrays
X=zeros(2,2,N+1); %true data
P=zeros(2,4,N+1); %filtration error covariance matrix
K=zeros(2,N);
X_f=zeros(2,N);

%Initial conditions
X(:,:,1)=[2,0;0,0]; 
P(:,1:2,1)=[10000 0; 0 10000];
Q=G*G.'*sigma2_a;
R=sigma2_n;

for i=2:N+1
    %Prediction of state vector at time i using i-1 measurements
    X(:,2,i-1)=Fi*X(:,1,i-1);
    X_f(:,i-1)=Fi^7*X(:,1,i-1);
    %Prediction error covariance matrix
    P(:,3:4,i-1)=Fi*P(:,1:2,i-1)*Fi.'+Q;
    %Filter gain, weight of residual
    K(:,i-1)=P(:,3:4,i-1)*H.'*(H*P(:,3:4,i-1)*H.'+R)^(-1);
    %Improved estimate by incorporating a new measurement
    X(:,1,i)=X(:,2,i-1)+K(:,i-1)*(Z(i-1)-H*X(:,2,i-1));
    %Filtration error covariance matrix
    P(:,1:2,i)=(eye(2)-K(:,i-1)*H)*P(:,3:4,i-1);
end

x_Kalman(1:N)=X(1,1,2:N+1);

%% 5) Building of plots of true trajectory, measurments, filtered esttimates of
%state vector
t=1:N; %array of steps
figure
plot(t,x,'g',t,Z,'r',t,x_Kalman,'k')
legend('true data','measurments','filtered estimates of state vector');
grid on
xlabel('Step')
ylabel('Data')
title('Applying of backward exponential mean')

figure
%Building of plot of changing filter gain and filtration error
subplot(2,1,1);
plot(t,K(1,:));
legend('K');
grid on
xlabel('Step')
ylabel('Filter gain K, filtration error')
title('changing of filter gain K');

%% 6)

P_sq(1:N)=sqrt(P(1,1,1:N));
subplot(2,1,2);
plot(t,P_sq);
legend('filtrating error');
grid on
xlabel('Step')
ylabel('Filtration error')
title('changing of filtration error');

%% 7-8-9)

%Estimation of dinamics of mean-sqaured error of estimation over
%observation interval
clear all

%Initial data
m=7; %number of steps ahead
M=500; %number of runs
N=200; %observation interval
sigma2_n=20^2; %variance of noise
sigma2_a=0.2^2; %variance of acceleraration
T=1; %period of step

%Creation of array of errors
Error_X=zeros(M,N); %array of errors of filtered estimate
Error_X_f=zeros(M,N); %array of errors of forecasts
Error_X_f7=zeros(M,N-6); %array of errors of forecasts
X_Kalman=zeros(6,N,M); %array of filtered data
x=zeros(2,N);

P(:,1:2,1)=[10000 0; 0 10000];

for i=1:M
    x=Generation_true_a(N,T,sigma2_n,sigma2_a);
    Z=x(2,:);
    X_Kalman(:,:,i)=Kalman_filter(Z,P,T,m,sigma2_n,sigma2_a); %Kalman filter
    Error_X(i,:)=(x(1,:)-X_Kalman(1,:,i)).^2; %errors of filtered estimate
    Error_X_f(i,:)=(x(1,:)-X_Kalman(2,:,i)).^2; %errors of forecasts 1-step
    Error_X_f7(i,:)=(x(1,7:N)-X_Kalman(3,1:N-6,i)).^2; %errors of forecasts 7-step
end

%Final average value of Error over M runs
Final_Error_X=sqrt(1/(M-1)*sum(Error_X));
Final_Error_X_f=sqrt(1/(M-1)*sum(Error_X_f));
Final_Error_X_f7=sqrt(1/(M-1)*sum(Error_X_f7));

%Building of plot of final errors of filtered estimate and errors of
%forecasts
t=1:N-2;
figure
%True error for filtration
subplot(2,1,1);
plot(t,Final_Error_X(3:N),t,X_Kalman(4,3:N,1));
legend('True estimation error','Filtration error covariance matrix');
grid on
xlabel('Step')
ylabel('Errors')
title('Final error of filtered estimates');

%Prediction error 1-step ahead
subplot(2,1,2);
plot(t,Final_Error_X_f(3:N),t,X_Kalman(5,3:N,1));
legend('True error of 1-step forecasts','Prediction error covariance matrix');
grid on
xlabel('Step')
ylabel('Errors')
title('Final error of 1-step forecasts');

%True error for filtration, prediction 1-step ahead and 7-steps ahead
t1=7:N;
figure
plot(t,Final_Error_X(3:N),t,Final_Error_X_f(3:N),t1,Final_Error_X_f7);
legend('True estimation error','True error of 1-step forecasts','True error of 7-step forecasts');
grid on
xlabel('Step')
ylabel('Errors')
title('Final error true estimates of 1-step and 7-step forecasts');

%Shoving how true estimation error, filtration error covariance matrix and
%filter gain approache to zero
t=1:N;
figure
subplot(2,1,1);
plot(t,Final_Error_X,t,X_Kalman(4,:,1));
legend('True estimation error','Filtration error covariance matrix');
grid on
xlabel('Step')
ylabel('Errors')
title('Errors changing');

subplot(2,1,2);
plot(t,X_Kalman(6,:,1));
legend('Filter gain');
grid on
xlabel('Step')
ylabel('K')
title('Filter gain changing');

%% 10) Make 𝑀=500  runs again, but with more accurate initial filtration
% error covariance matrix 

% NEW P(0,0)    
P_n(:,1:2,1)=[100 0; 0 100];

for i=1:M
    x_n=Generation_true_a(N,T,sigma2_n,sigma2_a);
    Z_n=x_n(2,:);
    X_Kalman_n(:,:,i)=Kalman_filter(Z_n,P_n,T,m,sigma2_n,sigma2_a); %Kalman filter
    Error_X_n(i,:)=(x(1,:)-X_Kalman_n(1,:,i)).^2; %errors of filtered estimate
    Error_X_f_n(i,:)=(x(1,:)-X_Kalman_n(2,:,i)).^2; %errors of forecasts 1-step
    Error_X_f7_n(i,:)=(x(1,7:N)-X_Kalman_n(3,1:N-6,i)).^2; %errors of forecasts 7-step
end

%Final average value of Error over M runs
Final_Error_X_n=sqrt(1/(M-1)*sum(Error_X_n));
Final_Error_X_f_n=sqrt(1/(M-1)*sum(Error_X_f_n));
Final_Error_X_f7_n=sqrt(1/(M-1)*sum(Error_X_f7_n));

clc
clear
close all

%% Part 1

% TRUE TRAJECTORY
steps=200; % size of trajectory
x=5*ones(steps,1); % initial coordinate
v=ones(steps,1); % initial velocity
t=1; % period
sigma_a=0.2; % variance of noise a 
sigma_eta=20; % variance of noise eta

for i=2:steps
    a=normrnd(0,1)*sigma_a; % generate random noise
    x(i)=x(i-1)+v(i-1)*t+a*t^2/2; % true trajectory of an object
    v(i)=v(i-1)+a*t; % calculate velocity
end

% MEASUREMENT OF COORDINATE X
for i=1:steps
    z(i)=x(i)+normrnd(0,1)*sigma_eta; % add noise to the trajectory
end
 
% SYSTEM AT STATE SPACE
F=[1 t;0 1]; % Transition matrix
G=[t^2/2;t]; % Input matrix
H=[1 0]; % Observation matrix

% KALMAN FILTER ALGORITHM DEVELOPMENT

X_ka=[2*ones(1,steps);zeros(1,steps)]; % state vector, describes full state of the system(coordinate&velocity)
Q=G*G'*sigma_a^2; % covariance matrix of state noise
P(1:2,1:2)=[10000 0;0 10000]; % Initial filtration error covariance matrix
R=sigma_eta^2; 
K=zeros(2,steps);
j=-1;

for i=2:steps
    X_ka(:,i)=F*X_ka(:,i-1); % construct kalman filter
    j=j+2;
    P_pred(1:2,j+2:j+3)=F*P(1:2,j:j+1)*F'+Q; % filling filtration error covariance matrix
    K(1:2,i-1)=P_pred(1:2,j+2:j+3)*H'/(H*P_pred(1:2,j+2:j+3)*H'+R); % FILTER GAIN
    X_ka(:,i)=X_ka(:,i)+K(1:2,i-1)*(z(i)-H*X_ka(:,i)); % Stationary Kalman filter
    P(1:2,j+2:j+3)=(eye(2)-K(1:2,i-1)*H)*P_pred(1:2,j+2:j+3); % error covarience matrix
end

% K(:,1)= [NaN,NaN]';

figure;
plot(1:steps,z,'r','Linewidth',1.2) % plot measurement
hold on
plot(1:steps,x,'b','Linewidth',1.2) % plot true data
plot(1:steps,X_ka(1,:),'k','Linewidth',1.2) % plot filtered data
grid on
title('Filtered data & True data & Measurements','Fontweight','bold');
xlabel('Steps','Fontweight','bold');
ylabel('Data','Fontweight','bold');
legend('Measurements Z_i','True data X_i','Filtered data','Fontweight','bold');

% GENERATE MEASUREMENTS WITH GAPS
for i=2:steps
 eps=normrnd(0,2);
 eps=sqrt(eps^2);
 if eps<=0.2
     z(i)=NaN;
 end 
end

figure
plot(1:steps,z,'r','Linewidth',1.2) % plot measurement
hold on
plot(1:steps,x,'b','Linewidth',1.2) % plot true data
grid on
title('Measurements with gaps','Fontweight','bold');
xlabel('Steps','Fontweight','bold');
ylabel('Data','Fontweight','bold');
legend('Measurements Z_i with gaps','True data X_i','Fontweight','bold');

% SYSTEM AT STATE SPACE
F=[1 t;0 1]; % Transition matrix
G=[t^2/2;t]; % Input matrix
H=[1 0]; % Observation matrix
% KALMAN FILTER ALGORITHM FOR TRAJECTORY WITH GAPS

X_ka=[2*ones(1,steps);zeros(1,steps)]; % state vector, describes full state of the system(coordinate&velocity)
Q=G*G'*sigma_a^2; % covariance matrix of state noise
P(1:2,1:2)=[10000 0;0 10000]; % Initial filtration error covariance matrix
R=sigma_eta^2; % covariance matrix of measurements noise  
K=zeros(2,steps);
j=-1;

% FIX NaN IN MEASUREMENTS & KALMAN FILTER DEVELOPMENT
 
for i=2:steps
 if isnan(z(i))
     z(i)=z(i-1);
 end 
 X_ka(:,i)=F*X_ka(:,i-1); % construct kalman filter
 j=j+2;
 P_pred(1:2,j+2:j+3)=F*P(1:2,j:j+1)*F'+Q; % filling iltration error covariance matrix
 K(1:2,i-1)=P_pred(1:2,j+2:j+3)*H'/(H*P_pred(1:2,j+2:j+3)*H'+R); % FILTER GAIN
 X_ka(:,i)=X_ka(:,i)+K(1:2,i-1)*(z(i)-H*X_ka(:,i)); % Stationary Kalman filter
 P(1:2,j+2:j+3)=(eye(2)-K(1:2,i-1)*H)*P_pred(1:2,j+2:j+3); % error covarience matrix
end

figure()
plot(1:steps,z,'Linewidth',1.2) % plot measurement
hold on
plot(1:steps,X_ka(1,:),'Linewidth',1.2) % plot filtered data
grid on;
title('Filtered data & Fixed measurements','Fontweight','bold');
xlabel('Steps','Fontweight','bold');
ylabel('Data','Fontweight','bold');
legend('Fixed measurements Z_i','Filtered data','Fontweight','bold');

%% ERRORS

prob = 0.3;
% prob = 0.5;
% prob = 0.7;

% 
% [final_error_x, final_error_pred, final_error_set] = err(prob, steps, x, v,...
%     t, sigma_a,sigma_eta);

for i=2:steps
 eps=normrnd(0,3);
 eps=sqrt(eps^2);
 if eps<=0.2
 z(i)=NaN;
 end 
end

for M=1:500
    for i=2:steps
        a=normrnd(0,1)*sigma_a;
        x(i)=x(i-1)+v(i-1)*t+a*t^2/2;
        v(i)=v(i-1)+a*t;
    end
    z=x;
    for i=1:steps
        zeta=normrnd(0,3);
        zeta=sqrt(zeta^2);
        if zeta>prob
            z(i)=x(i)+normrnd(0,1)*sigma_eta;
        else
            z(i)=NaN;
        end
    end
    for i=2:steps
        X_ka(:,i)=F*X_ka(:,i-1);
        X_set(:,i)=F^7*X_ka(:,i);
        X_pred(:,i)=F*X_ka(:,i-1);
        predict(1:2,i*2-1:i*2)=F*P(1:2,i*2-3:i*2-2)*F'+Q;
        k(1:2,i-1)=predict(1:2,i*2-1:i*2)*H'/(H*predict(1:2,i*2-1:i*2)*H'+R);
        if isnan(z(i))
            X_ka(:,i)=X_ka(:,i);
            P(1:2,i*2-1:i*2)=P(1:2,i*2-1:i*2);
        else
            X_ka(:,i)=X_ka(:,i)+k(1:2,i-1)*(z(i)-H*X_ka(:,i));
            P(1:2,i*2-1:i*2)=(eye(2)-k(1:2,i-1)*H)*predict(1:2,i*2-1:i*2);
        end
    end
    error_x(:,M)=(x(3:end)-X_ka(1,3:end)').^2;
    error_pred(:,M)=(x(3:end)-X_pred(1,3:end)').^2;
    error_set(:,M)=(x(7:end)-X_set(1,2:end-5)').^2;
end

final_error_x_03=sqrt(1/(M-1)*sum(error_x,2));
final_error_pred=sqrt(1/(M-1)*sum(error_pred,2));
final_error_set=sqrt(1/(M-1)*sum(error_set,2));

figure()
plot(1:steps-6,final_error_x_03(1:end-4),'m','Linewidth',1.2)
hold on
plot(1:steps-6,final_error_pred(1:end-4),'c','Linewidth',1.2)
plot(1:steps-6,final_error_set,'k','Linewidth',1.2)
grid on; grid minor
xlabel('Step', 'FontSize', 30)
ylabel('Errors', 'FontSize', 30)
legend('True estimation error','True error of 1-step prediction with P=0.5', ...
    'True error of 7-step prediction with P=0.5', 'FontSize', 30)

%% P = 0.5

prob = 0.5;

for i=2:steps
 eps=normrnd(0,3);
 eps=sqrt(eps^2);
 if eps<=0.2
 z(i)=NaN;
 end 
end

for M=1:500
    for i=2:steps
        a=normrnd(0,1)*sigma_a;
        x(i)=x(i-1)+v(i-1)*t+a*t^2/2;
        v(i)=v(i-1)+a*t;
    end
    z=x;
    for i=1:steps
        zeta=normrnd(0,3);
        zeta=sqrt(zeta^2);
        if zeta>prob
            z(i)=x(i)+normrnd(0,1)*sigma_eta;
        else
            z(i)=NaN;
        end
    end
    for i=2:steps
        X_ka(:,i)=F*X_ka(:,i-1);
        X_set(:,i)=F^7*X_ka(:,i);
        X_pred(:,i)=F*X_ka(:,i-1);
        predict(1:2,i*2-1:i*2)=F*P(1:2,i*2-3:i*2-2)*F'+Q;
        k(1:2,i-1)=predict(1:2,i*2-1:i*2)*H'/(H*predict(1:2,i*2-1:i*2)*H'+R);
        if isnan(z(i))
            X_ka(:,i)=X_ka(:,i);
            P(1:2,i*2-1:i*2)=P(1:2,i*2-1:i*2);
        else
            X_ka(:,i)=X_ka(:,i)+k(1:2,i-1)*(z(i)-H*X_ka(:,i));
            P(1:2,i*2-1:i*2)=(eye(2)-k(1:2,i-1)*H)*predict(1:2,i*2-1:i*2);
        end
    end
    error_x(:,M)=(x(3:end)-X_ka(1,3:end)').^2;
    error_pred(:,M)=(x(3:end)-X_pred(1,3:end)').^2;
    error_set(:,M)=(x(7:end)-X_set(1,2:end-5)').^2;
end

final_error_x_05=sqrt(1/(M-1)*sum(error_x,2));

%% P = 0.7

prob = 0.7;

for M=1:500
    for i=2:steps
        a=normrnd(0,1)*sigma_a;
        x(i)=x(i-1)+v(i-1)*t+a*t^2/2;
        v(i)=v(i-1)+a*t;
    end
    z=x;
    for i=1:steps
        zeta=normrnd(0,3);
        zeta=sqrt(zeta^2);
        if zeta>prob
            z(i)=x(i)+normrnd(0,1)*sigma_eta;
        else
            z(i)=NaN;
        end
    end
    for i=2:steps
        X_ka(:,i)=F*X_ka(:,i-1);
        X_set(:,i)=F^7*X_ka(:,i);
        X_pred(:,i)=F*X_ka(:,i-1);
        predict(1:2,i*2-1:i*2)=F*P(1:2,i*2-3:i*2-2)*F'+Q;
        k(1:2,i-1)=predict(1:2,i*2-1:i*2)*H'/(H*predict(1:2,i*2-1:i*2)*H'+R);
        if isnan(z(i))
            X_ka(:,i)=X_ka(:,i);
            P(1:2,i*2-1:i*2)=P(1:2,i*2-1:i*2);
        else
            X_ka(:,i)=X_ka(:,i)+k(1:2,i-1)*(z(i)-H*X_ka(:,i));
            P(1:2,i*2-1:i*2)=(eye(2)-k(1:2,i-1)*H)*predict(1:2,i*2-1:i*2);
        end
    end
    error_x(:,M)=(x(3:end)-X_ka(1,3:end)').^2;
    error_pred(:,M)=(x(3:end)-X_pred(1,3:end)').^2;
    error_set(:,M)=(x(7:end)-X_set(1,2:end-5)').^2;
end

final_error_x_07=sqrt(1/(M-1)*sum(error_x,2));

%% Final plot

figure()
plot(1:steps-6,final_error_x_03(1:end-4),'m','Linewidth',1.2)
hold on
plot(1:steps-6,final_error_x_05(1:end-4),'c','Linewidth',1.2)
plot(1:steps-6,final_error_x_07(1:end-4),'k','Linewidth',1.2)
grid on; grid minor
xlabel('Step', 'FontSize', 30)
ylabel('Errors', 'FontSize', 30)
legend('True estimation error with P=0.3','True estimation error with P=0.5', ...
    'True estimation error with P=0.7', 'FontSize', 30)
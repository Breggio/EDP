function [final_error_x, final_error_pred, final_error_set] = err(prob, steps, x,...
    v, t, sigma_a,sigma_eta)

F=[1 t;0 1]; % Transition matrix
G=[t^2/2;t]; % Input matrix
H=[1 0]; % Observation matrix

X_ka=[2*ones(1,steps);zeros(1,steps)]; % state vector, describes full state of the system(coordinate&velocity)
Q=G*G'*sigma_a^2; % covariance matrix of state noise
P(1:2,1:2)=[10000 0;0 10000]; % Initial filtration error covariance matrix
R=sigma_eta^2; % covariance matrix of measurements noise  
K=zeros(2,steps);

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

final_error_x=sqrt(1/(M-1)*sum(error_x,2));
final_error_pred=sqrt(1/(M-1)*sum(error_pred,2));
final_error_set=sqrt(1/(M-1)*sum(error_set,2));
end


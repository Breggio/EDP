function [x, z, x_Kalman] = tracking(N, T, sigma2_a, sigma2_n, Prob)

eta = randn*sqrt(sigma2_n); % Random noise of measurments

x(1) = 5;
V(1) = 1;
z(1) = x(1) + eta; %the first measurment

%State matrixes
phi = [1 T; 0 1];
G  = [T/2; T];    % input matrix
H  = [1 0];       % observation matrix

% Kalman filter 
X = zeros(2, 2, N + 1); %true data
P = zeros(2, 4, N + 1); %filtration error covariance matrix
K = zeros(2, N);
X_f = zeros(2, N);

% Initial conditions
X(:,:,1) = [2, 0; 0, 0];
P(:,1:2,1) = [10000 0; 0 10000];
Q = G*G.'*sigma2_a; 
R = sigma2_n;
m = 7; % steps ahead

for i=2:N+1

    a = randn*sqrt(sigma2_a); %normally distributed random acceleration
    eta = randn*sqrt(sigma2_n); %random noise of measurments
    V(i) = V(i - 1) + a*T;
    x(i) = x(i - 1) + V(i - 1)*T + a*T^2/2;

    csi = rand;

    if csi <= Prob
        z(i) = NaN;
        X(:,:,i) = X(:,:,i-1);
        P(:,:,i) = P(:,:,i-1);
        X_f(:,i) = phi^7*X(:,1,i);
    
    elseif csi > Prob
        z(i) = x(i) + eta;
        X_f(:,i-1) = phi^7*X(:,1,i-1);
        X(:,2,i-1) = phi*X(:,1,i-1);
        P(:,3:4,i-1) = phi*P(:,1:2,i-1)*phi.'+Q;
        K(:,i-1) = P(:,3:4,i-1)*H.'*(H*P(:,3:4,i-1)*H.'+R) ^(-1);

        if isnan(z(i-1)) == 1 
            counter=i ;
            znew=z(i-1);
            
            while isnan(znew) == 1 
                znew=z ( counter-1); 
                counter=counter-1;
            end
        elseif isnan(z(i-1)) == 0 
            znew=z(i-1);
        end
        
        X(:,1,i) = X(:,2,i-1)+K(:,i-1)*(znew-H*X(:,2, i-1) ) ;
        P(:,1:2,i) = (eye(2)-K(:,i-1)*H)*P(:,3:4,i-1);

    end
end

x_Kalman(1:N) = X(1,1,2:N+1);

end


function [X, P, D_K, b_K, condition_num, K] = Kalman_filter_determ(Z_c, Z_p, T, sigma_D, sigma_beta)
    N = length(Z_c);
    Fi = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1]; %transition matrix
    H = [1 0 0 0; 0 0 1 0]; %observation matrix

    %Kalman filter development
    %Creating of arrays
    X = zeros(4,2,N+1);   %filtered and smoothed data
    P = zeros(4,8,N+1);   %filtration and prediction error covariance matrix
    K = zeros(4,2,N);     %filter gain
    D_K = zeros(2,N);     %array of filtered and extrapolated range
    b_K = zeros(2,N);     %array of filtered and extrapolated azimuth
    lambda1 = zeros(1,N); %array of first eigenvalues
    lambda2 = zeros(1,N); %array of second eigenvalues
    condition_num = zeros(1,N); %array of condition numbers

    %Initial conditions
    X(:,:,1)=[40000,0; -20,0; 40000,0; -20,0]; 
    P(:,1:4,1)=[10^10 0 0 0; 0 10^10 0 0; 0 0 10^10 0; 0 0 0 10^10];
    R = zeros(2, 2);
    for i=2:N+1 %0 1 2 3 .... N N+1
        R(1, 1) = sin(Z_p(2,i - 1))^2*sigma_D^2 + (Z_p(1,i - 1))^2*cos(Z_p(2,i - 1))^2*sigma_beta^2;
        R(2, 2) = cos(Z_p(2,i - 1))^2*sigma_D^2 + (Z_p(1,i - 1))^2*sin(Z_p(2,i - 1))^2*sigma_beta^2;
        R(1, 2) = sin(Z_p(2,i - 1))*cos(Z_p(2,i - 1))*(sigma_D^2 - (Z_p(1,i - 1))^2*sigma_beta^2);
        R(2, 1) = sin(Z_p(2,i - 1))*cos(Z_p(2,i - 1))*(sigma_D^2 - (Z_p(1,i - 1))^2*sigma_beta^2);
        %Prediction of state vector at time i using i-1 measurements
        X(:,2,i - 1) = Fi*X(:,1,i-1);
        %Prediction error covariance matrix
        P(:,5:8,i - 1) = Fi*P(:,1:4,i-1)*Fi.';
        %Filter gain, weight of residual
        K(:,:,i-1) = P(:,5:8,i-1)*H.'*(H*P(:,5:8,i-1)*H.'+R)^(-1);
        %Improved estimate by incorporating a new measurement
        X(:,1,i) = X(:,2,i-1) + K(:,:,i-1)*(Z_c(i-1) - H*X(:,2,i-1));
        %Filtration error covariance matrix
        P(:,1:4,i) = (eye(4) - K(:,:,i-1)*H)*P(:,5:8,i-1);
        %Evaluation of predicted and extrapolated range and azimuth
        D_K(1,i-1) = sqrt(X(1,1,i)^2 + X(3,1,i)^2);
        D_K(2,i-1) = sqrt(X(1,2,i-1)^2 + X(3,2,i-1)^2);
        b_K(1,i-1) = atan(X(1,1,i)/X(3,1,i));
        b_K(2,i-1) = atan(X(1,2,i-1)/X(3,2,i-1));
        lambda1(i-1) = sigma_D^2; %const 400
        lambda2(i-1) = Z_p(1,i-1)^2*sigma_beta^2; %descrise 10^4
        if lambda1(i-1)>lambda2(i-1)
            condition_num(i-1) = lambda1(i-1)/lambda2(i-1);
        else
            condition_num(i-1) = lambda2(i-1)/lambda1(i-1);
        end
    end
end
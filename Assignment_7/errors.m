function [Final_Error_X, Final_Error_X_f, Final_Error_X_f7] = errors(M, N,sigma2_n, sigma2_a, T, Prob)

% Initialization
Error_X    = zeros(M, N);      % array of errors of filtered estimate
Error_X_f  = zeros(M, N);      % array of errors of forecasts
Error_X_f7 = zeros(M, N - 6);  % array of errors of forecasts

for i = 1:M
    [x, ~, X_Kalman(:,:,i)] = tracking(N, T, sigma2_a, sigma2_n, Prob);
    Error_X(i,:)    = (x(1,:) - X_Kalman(1,:,i)).^2; % errors of filtered estimate
    Error_X_f(i,:)  = (x(1,:) - X_Kalman(2,:,i)).^2; % errors of forecasts 1-step
    Error_X_f7(i,:) = (x(1,7:N) - X_Kalman(3,1:N-6,i)).^2; % errors of forecasts 7-step
end

%Final average value of Error over M runs
Final_Error_X    = sqrt(1/(M-1)*sum(Error_X));
Final_Error_X_f  = sqrt(1/(M-1)*sum(Error_X_f));
Final_Error_X_f7 = sqrt(1/(M-1)*sum(Error_X_f7));

end
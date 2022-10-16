function [Z_sm, P_sm] = kalman_sm(data, P_pred, P, T)
    N = length(data);
    Fi = [1 T; 0 1];
    A = zeros(2,2,N - 1);  %coefficient
	P_sm = zeros(2, 2, N); %smoothing error covariance matrix
	Z_sm = zeros(2,N);     %smoothed data
    
    P_sm(:,:,N) = P(:,:,N);
	Z_sm(:,N) = data(:,N);   
    for i = N - 1:-1:1
        A(:,:,i) = P(:,:,i)*Fi.'/P_pred(:,:,i);
        P_sm(:,:,i) = P(:,:,i) + A(:,:,i)*(P_sm(:,:,i+1) - P_pred(:,:,i))*A(:,:,i).';
        Z_sm(:,i) = data(:,i) + A(:,:,i)*(Z_sm(:,i+1) - Fi*data(:,i));
    end
end


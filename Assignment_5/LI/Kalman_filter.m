function X_Kalman = Kalman_filter(Z,T,m,sigma2_n,sigma2_a)
N=length(Z);
Fi=[1 T; 0 1]; %transition matrix
G=[T/2;T];     %input matrix
H=[1 0];       %observation matrix

%Kalman filter development
%Creating of arrays
X=zeros(2,2,N+1); %true data
P=zeros(2,4,N+1); %filtration error covariance matrix
K=zeros(2,N);
X_f=zeros(2,N);

%Initial conditions
X(:,:,1) = [2,0;0,0]; 
P(:,1:2,1) = [10000 0; 0 10000];
Q = G*G.'*sigma2_a;
R = sigma2_n;

for i=2:N+1
    %Prediction of state vector at time i using i-1 measurements
    X(:,2,i-1)=Fi*X(:,1,i-1);
    X_f(:,i-1)=Fi^m*X(:,1,i-1);
    %Prediction error covariance matrix
    P(:,3:4,i-1)=Fi*P(:,1:2,i-1)*Fi.'+Q;
    %Filter gain, weight of residual
    K(:,i-1)=P(:,3:4,i-1)*H.'*(H*P(:,3:4,i-1)*H.'+R)^(-1);
    %Improved estimate by incorporating a new measurement
    X(:,1,i)=X(:,2,i-1)+K(:,i-1)*(Z(i-1)-H*X(:,2,i-1));
    %Filtration error covariance matrix
    P(:,1:2,i)=(eye(2)-K(:,i-1)*H)*P(:,3:4,i-1);
end

X_Kalman(1,1:N)=X(1,1,2:N+1); %filtered data
X_Kalman(2,1:N)=X(1,2,1:N); %predictions
X_Kalman(3,1:N-6)=X_f(1,1:N-6); %predictions 7 steps ahead
X_Kalman(4,1:N)=sqrt(P(1,1,2:N+1)); %error of filtration
X_Kalman(5,1:N)=sqrt(P(1,3,1:N)); %error of prediction
X_Kalman(6,1:N)=K(1,:); %filter gain
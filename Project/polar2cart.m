function [Z_c_1, Z_c_2, Z_c_3] = polar2cart(N, z1, z2, z3, sigma_d, sigma_beta)

d = [z1(:,1), z2(:,1), z3(:,1)];
beta = [z1(:,2), z2(:,2), z3(:,2)];

%true_polar = [d; beta];

% Initialization

d_m_1 = zeros(1,N);
d_m_2 = zeros(1,N);
d_m_3 = zeros(1,N);
b_m_1 = zeros(1,N);
b_m_2 = zeros(1,N);
b_m_3 = zeros(1,N);
% 
% x_m_1 = zeros(1,N);
% x_m_2 = zeros(1,N);
% x_m_3 = zeros(1,N);
% 
% y_m_1 = zeros(1,N);
% y_m_2 = zeros(1,N);
% y_m_3 = zeros(1,N);
% 
% Z_c_1 = zeros(1,N);
% Z_c_2 = zeros(1,N);
% Z_c_3 = zeros(1,N);
% 
% Z_p_1 = zeros(1,N);
% Z_p_2 = zeros(1,N);
% Z_p_3 = zeros(1,N);

for i = 1:N
    a=randn;
    d_m_1(i) = d(1,i) + a*sigma_d;
    b_m_1(i) = beta(1,i) + a*sigma_beta;
    d_m_2(i) = d(2,i) +a*sigma_d;
    b_m_2(i) = beta(2,i) + a*sigma_beta;
    d_m_3(i) = d(3,i) + a*sigma_d;
    b_m_3(i) = beta(3,i) + a*sigma_beta;

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


end


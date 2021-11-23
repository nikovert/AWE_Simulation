function vec_A = transformFromBtoA( alpha, beta,  vec_B )


c_a = cos(alpha);
s_a = sin(alpha);

c_b = cos(beta);
s_b = sin(beta);

M_B2A = zeros(3,3);

M_B2A(1,1) = c_a*c_b;
M_B2A(2,1) = -c_a*s_b;
M_B2A(3,1) = -s_a;
M_B2A(1,2) = s_b;
M_B2A(2,2) = c_b;
M_B2A(3,2) = 0;
M_B2A(1,3) = s_a*c_b;
M_B2A(2,3) = -s_a*s_b;
M_B2A(3,3) = c_a;

vec_A = M_B2A * vec_B; 
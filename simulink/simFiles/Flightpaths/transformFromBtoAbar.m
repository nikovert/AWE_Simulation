function vec_Abar = transformFromBtoAbar( mu_a, alpha, beta,  vec_B )

c_m = cos(mu_a);
s_m = sin(mu_a);

c_a = cos(alpha);
s_a = sin(alpha);

c_b = cos(beta);
s_b = sin(beta);

M_B2ABar = zeros(3,3);

M_B2ABar(1,1) = c_a*c_b;
M_B2ABar(1,2) = s_b;
M_B2ABar(1,3) = s_a*c_b;
M_B2ABar(2,1) = -c_a*s_b*c_m+s_a*s_m;
M_B2ABar(2,2) = c_b*c_m;
M_B2ABar(2,3) = -s_a*s_b*c_m-s_m*c_a;
M_B2ABar(3,1) = -c_a*s_b*s_m-s_a*c_m;
M_B2ABar(3,2) = c_b*s_m;
M_B2ABar(3,3) = -s_a*s_b*s_m+c_a*c_m;

vec_Abar = M_B2ABar * vec_B; 
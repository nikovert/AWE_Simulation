function vec_B = transformFromAbar2B( mu_a, alpha, beta,  vec_Abar )

c_m = cos(mu_a);
s_m = sin(mu_a);

c_a = cos(alpha);
s_a = sin(alpha);

c_b = cos(beta);
s_b = sin(beta);

M_ABar2B = zeros(3,3);

M_ABar2B(1,1) = c_a*c_b;
M_ABar2B(2,1) = s_b;
M_ABar2B(3,1) = s_a*c_b;
M_ABar2B(1,2) = -c_a*s_b*c_m+s_a*s_m;
M_ABar2B(2,2) = c_b*c_m;
M_ABar2B(3,2) = -s_a*s_b*c_m-s_m*c_a;
M_ABar2B(1,3) = -c_a*s_b*s_m-s_a*c_m;
M_ABar2B(2,3) = c_b*s_m;
M_ABar2B(3,3) = -s_a*s_b*s_m+c_a*c_m;

vec_B = M_ABar2B * vec_Abar; 
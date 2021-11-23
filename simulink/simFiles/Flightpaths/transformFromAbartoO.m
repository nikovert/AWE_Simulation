function vec_O = transformFromAbartoO( chi_a,gamma_a,  vec_Abar )

c_chia = cos(chi_a);
s_chia = sin(chi_a);

c_gama = cos(gamma_a);
s_gama = sin(gamma_a);

M_ABar2O = zeros(3,3);

M_ABar2O(1,1) = c_chia*c_gama;
M_ABar2O(2,1) = s_chia*c_gama;
M_ABar2O(3,1) = -s_gama;
M_ABar2O(1,2) = -s_chia;
M_ABar2O(2,2) = c_chia;
M_ABar2O(3,2) = 0;
M_ABar2O(1,3) = c_chia*s_gama;
M_ABar2O(2,3) = s_chia*s_gama;
M_ABar2O(3,3) = c_gama;

vec_O = M_ABar2O * vec_Abar; 
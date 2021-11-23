function vec_Abar = transformFromO2Abar( chi_a, gamma_a, vec_O )

M_AbarO = [cos(chi_a)*cos(gamma_a), sin(chi_a)*cos(gamma_a), -sin(gamma_a); 
    -sin(chi_a), cos(chi_a), 0; 
    cos(chi_a)*sin(gamma_a), sin(chi_a)*sin(gamma_a), cos(gamma_a)];
vec_Abar = M_AbarO*vec_O; 
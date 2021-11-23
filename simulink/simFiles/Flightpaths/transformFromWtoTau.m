function vec_tau = transformFromWtoTau(  lamb,phi, vec_W )

M_tauW = [-sin(phi)*cos(lamb), -sin(phi)*sin(lamb), cos(phi);
    -sin(lamb), cos(lamb), 0;
    -cos(phi)*cos(lamb), -cos(phi)*sin(lamb), -sin(phi)];

vec_tau = M_tauW * vec_W; 
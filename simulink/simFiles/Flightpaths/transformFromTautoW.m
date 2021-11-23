function vec_W = transformFromTautoW( lamb,phi,  vec_tau )

M_Wtau = [-sin(phi)*cos(lamb),-sin(lamb),-cos(phi)*cos(lamb);
    -sin(phi)*sin(lamb), cos(lamb), -cos(phi)*sin(lamb);
    cos(phi), 0, -sin(phi)];


vec_W = M_Wtau * vec_tau; 
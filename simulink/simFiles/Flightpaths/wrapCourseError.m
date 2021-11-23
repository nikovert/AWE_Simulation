function e_course = wrapCourseError(chi_ref, chi)
%#codegen
% The vectors are lying in the tangential plane attached to the current
% position on the small earth.
p_bearing_vec =   [cos(chi);         sin(chi);         0];
ref_bearing_vec = [cos(chi_ref); sin(chi_ref); 0];

e_z = cross( ref_bearing_vec, p_bearing_vec ); 

dot_product = dot( ref_bearing_vec, p_bearing_vec ); 
if dot_product > 1 
    dot_product = 1; 
end
if dot_product < -1 
    dot_product = -1; 
end
e_course = - sign( e_z(3) ) * acos( dot_product );

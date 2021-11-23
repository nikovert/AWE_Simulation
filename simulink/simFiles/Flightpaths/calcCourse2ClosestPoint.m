function [ chi_K_des, bearing_vec_W, bearing_sign ] = calcCourse2ClosestPoint( p_VT_W, p_kite_W )
% This function return the desired course in the tangential plane.
%================== Function Description ================================
% Version: 1.0 
% Author: Sebastian Rapp 
% Last Update: 28.08.2017
% Inputs: position of the kite and the virtual target 
% Output: the desired course and the vector pointing to the desired point
%========================================================================
bearing_sign = 0;
p_kite_W = p_kite_W/norm(p_kite_W);
p_VT_W = p_VT_W/norm(p_VT_W); 
e_1 = [p_VT_W(2)*p_kite_W(3)-p_VT_W(3)*p_kite_W(2); 
       p_VT_W(3)*p_kite_W(1)-p_VT_W(1)*p_kite_W(3); 
       p_VT_W(1)*p_kite_W(2)-p_VT_W(2)*p_kite_W(1)];
if norm(e_1) < 1e-9
    chi_K_des = 0;
    bearing_vec_W = [0;0;0];
else

e_1 = e_1/norm(e_1);

% Desired direction in world frame coordinates
bearing_vec_W = [p_kite_W(2)*e_1(3)-p_kite_W(3)*e_1(2); 
       p_kite_W(3)*e_1(1)-p_kite_W(1)*e_1(3); 
       p_kite_W(1)*e_1(2)-p_kite_W(2)*e_1(1)];
bearing_vec_W = bearing_vec_W/norm(bearing_vec_W); 

% Extract latitude/longitude with respect to the ground station
lat_kite = asin( p_kite_W(3) / norm(p_kite_W) ); 
long_kite = atan2( p_kite_W(2), p_kite_W(1) ); 

% Orthonormal basis of the tether frame
e_xT_W = [-sin(lat_kite)*cos(long_kite),-sin(lat_kite)*sin(long_kite), cos(lat_kite)]; 
e_zT_W = -p_kite_W/norm(p_kite_W); 
e_yT_W = [e_zT_W(2)*e_xT_W(3)-e_zT_W(3)*e_xT_W(2),e_zT_W(3)*e_xT_W(1)-e_zT_W(1)*e_xT_W(3),e_zT_W(1)*e_xT_W(2)-e_zT_W(2)*e_xT_W(1)];

% Calculate course
bearing_sign = sign( e_yT_W * bearing_vec_W ); % bearing_vec_T_y
chi_K_des =  bearing_sign * acos(  max( min( e_xT_W * bearing_vec_W, 1), -1 ) );% third component of cross(e_bearing, e_xT) determines sign
 
end
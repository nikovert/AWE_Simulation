function vrot = doRodriguesRotation(pos_W, p_target_W, v)

theta = acos( max( min( pos_W'*p_target_W / norm(pos_W) / norm( p_target_W ), 1),-1) );
k = cross( p_target_W, pos_W );

if abs(theta) < 1e-12 || norm(k) < 1e-12
    vrot = v;
else
    k = k/norm(k);
    vrot = v*cos(theta)+cross(k,v)*sin(theta)+k * (k'*v)*(1-cos(theta));
end
end
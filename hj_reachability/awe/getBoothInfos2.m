function [t,DtDs,L, dLds] = getBoothInfos2(s_old,Lbooth, direction)
    %GETBOOTHINFOS2 Calculate the tangent at the solution and its derivative
    % -1: down in the middle
    a = Lbooth.a; 
    b = Lbooth.b;
    L = { b .* sin(s_old) ./( 1+(a./b .* cos(s_old)).^2 );
          a .* sin(s_old).*cos(s_old) ./ ( 1+(a./b .* cos(s_old)).^2 ) } ;

    dLds = { ( b.^3 .* cos(s_old).*(2*a.^2-a.^2 .* cos(s_old).^2+b.^2)./(a.^2 .* cos(s_old).^2+b.^2).^2 );
        ( (cos(s_old).^2 .* (a.^3 .* b.^2+2*a .* b.^4) - a .* b.^4)./(a.^2 .* cos(s_old).^2+b.^2).^2 )};  
 
    d2Lds2 =  { ( -( (a.^4.*b.^3.*sin(s_old).^5-b.^3.*sin(s_old).*(5*a.^4+4*(a.*b).^2-b.^4)+b.^3.*sin(s_old).^3 .*(4*a.^4+6*(a.*b).^2 ))./(b.^2-a.^2.*(sin(s_old).^2 - 1 ) ).^3 ) );
                ( (2*a.*b.^2.*cos(s_old).^3.*sin(s_old).*(a.^4+2*(a.*b).^2)-a.*b.^2.*cos(s_old).*sin(s_old).*(3*(a.*b).^2+2*b.^4)*2)./(a.^2.*cos(s_old).^2 + b.^2 ).^3 )};


    s_lambda = sin( L{1}  ); 
    s_phi = sin( L{2} );
    c_lambda = cos( L{1} ); 
    c_phi = cos( L{2} ); 

    % q = {c_lambda.*c_phi; s_lambda.*c_phi; s_phi};

    dqdlambda = {-s_lambda.*c_phi; c_lambda.*c_phi; zeros(size(s_phi))};
    dqdphi = {-c_lambda.*s_phi; -s_lambda.*s_phi; c_phi};

    t = {direction*(dqdlambda{1} .* dLds{1} + dqdphi{1} .* dLds{2}); ...
         direction*(dqdlambda{2} .* dLds{1} + dqdphi{2} .* dLds{2}); ...
         direction*(dqdlambda{3} .* dLds{1} + dqdphi{3} .* dLds{2})};

    dtdlambda = {-c_lambda .* c_phi .* dLds{1} + s_lambda .* s_phi .* dLds{2}; ...
                 -s_lambda .* c_phi .* dLds{1} - s_phi .* c_lambda .* dLds{2}; ...
                0};

    dtdphi = {s_lambda .* s_phi .* dLds{1} - c_phi .* c_lambda .* dLds{2}; 
             -c_lambda .* s_phi .* dLds{1} - c_phi .* s_lambda .* dLds{2}; 
             -s_phi .* dLds{2}};

    dtds = {dqdlambda{1} .* d2Lds2{1} + dqdphi{1} .* d2Lds2{2}; ...
            dqdlambda{2} .* d2Lds2{1} + dqdphi{2} .* d2Lds2{2}; ...
            dqdlambda{3} .* d2Lds2{1} + dqdphi{3} .* d2Lds2{2}};
    
    % The negative sign cancels out 
    DtDs = {dtdlambda{1} .* dLds{1} + dtdphi{1} .* dLds{2} + dtds{1}; ...
            dtdlambda{2} .* dLds{1} + dtdphi{2} .* dLds{2} + dtds{2}; ...
                                      dtdphi{3} .* dLds{2} + dtds{3}};
end
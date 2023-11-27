function [array_fac] = array_factor(theta_vec, d, lambda,M)


psi = 2*pi/lambda*d*sind(theta_vec);
array_fac = exp(-1j*(M-1)./(2*psi)).*sin(M*psi/2)./(sin(psi/2));
    
return
